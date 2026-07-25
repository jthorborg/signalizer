#!/usr/bin/env bash
#
# Disassemble every code symbol whose (mangled) name contains a given
# substring, one .asm file per symbol, for tracking codegen changes
# through git over time.
#
# Usage: ./gold_disasm.sh <symbol-substring> [path-to-binary]
#
# For each matching symbol, writes:
#   gold/<arch>/<symbol-substring>/<mangled-name>.asm
# with the c++filt-demangled name as the first line, followed by the
# disassembly of that symbol.
#
# Note: nm matches are filtered to code symbols (T/t/W/w/I/i) since
# --disassemble-symbols only produces output for those; data symbols
# (e.g. local statics) that happen to share the substring are skipped.
#
# Note: unlike the Linux script (gold_disasm.sh in ../LinuxMakefile, which
# uses GNU binutils), this targets Apple's LLVM-based objdump/nm, which
# differ in a few important ways:
#   - the flag is --disassemble-symbols=, not --disassemble=
#   - -M intel doesn't select Intel syntax here (-M is aliased to
#     --disassembler-options); Intel syntax is --x86-asm-syntax=intel,
#     and only applies to x86 slices (arm64 has no such switch).
#   - --disassemble-symbols= takes a *comma-separated list*. Templated
#     C++ symbols routinely demangle to signatures containing commas
#     (e.g. "CComplexResonator<float, 2ul>"), which silently splits the
#     match into garbage fragments and matches nothing. There's no
#     escaping that helps. So, unlike the Linux script, we deliberately
#     do NOT pass -C/--demangle to objdump for the match: we match on the
#     raw mangled name (which contains no commas or shell metacharacters)
#     and only use c++filt to produce the human-readable header comment.
#     This means callee symbols referenced in the disassembly body stay
#     mangled, which is a worthwhile trade for actually working.
#
# Note: only Debug builds are usable here. Release builds LTO the archive
# members down to LLVM bitcode (no machine code to disassemble until final
# link), and the linked Release binary is stripped of the symbols we'd
# need to match on. Debug is unstripped and has real machine code, so it's
# the default target.
#
# Universal (multi-arch) binaries are thinned per-slice with lipo, and
# each arch's symbols/output are kept in their own gold/<arch>/ directory,
# since codegen differs by architecture.

set -euo pipefail

usage() {
    echo "Usage: $0 <symbol-substring> [path-to-binary]" >&2
    exit 1
}

[[ $# -ge 1 && $# -le 2 ]] || usage

substring=$1
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
binary=${2:-"$script_dir/build/Debug/Signalizer.app/Contents/MacOS/Signalizer"}

[[ -f "$binary" ]] || { echo "error: binary not found: $binary" >&2; exit 1; }

# Keep filenames well under typical filesystem limits (255 bytes); some
# heavily-templated symbols (e.g. JUCE's X11 symbol loader) mangle to
# well over a thousand characters.
max_name_len=200

archs=$(lipo -archs "$binary")

tmp_dir=""
cleanup() { [[ -n "$tmp_dir" ]] && rm -rf "$tmp_dir"; }
trap cleanup EXIT

total_written=0

for arch in $archs; do
    slice="$binary"
    if [[ $(echo "$archs" | wc -w) -gt 1 ]]; then
        tmp_dir=${tmp_dir:-$(mktemp -d)}
        slice="$tmp_dir/$arch"
        lipo -thin "$arch" -output "$slice" "$binary"
    fi

    out_dir="$script_dir/gold/$arch/$substring"

    symbols=()
    while IFS= read -r sym; do
        symbols+=("$sym")
    done < <(
        nm "$slice" | awk '$2 ~ /^[TtWwIi]$/ { print $3 }' | grep -F -- "$substring" | sort -u
    )

    if [[ ${#symbols[@]} -eq 0 ]]; then
        echo "no code symbols matching '$substring' found in $binary (arch $arch)" >&2
        continue
    fi

    mkdir -p "$out_dir"

    asm_flags=(-d -S -l --no-show-raw-insn)
    if [[ $arch == x86_64 || $arch == i386 ]]; then
        asm_flags+=(--x86-asm-syntax=intel)
    fi

    for sym in "${symbols[@]}"; do
        if (( ${#sym} > max_name_len )); then
            if command -v sha1sum >/dev/null 2>&1; then
                hash=$(printf '%s' "$sym" | sha1sum | cut -c1-10)
            else
                hash=$(printf '%s' "$sym" | shasum -a 1 | cut -c1-10)
            fi
            fname="${sym:0:$((max_name_len - 11))}_$hash.asm"
        else
            fname="$sym.asm"
        fi
        out_file="$out_dir/$fname"

        demangled=$(c++filt "$sym")

        echo "disassembling $sym -> gold/$arch/$substring/$fname"
        {
            echo "; $demangled"
            objdump "${asm_flags[@]}" --disassemble-symbols="$sym" "$slice"
        } > "$out_file"
    done

    echo "wrote ${#symbols[@]} file(s) to $out_dir"
    total_written=$((total_written + ${#symbols[@]}))
done

if [[ $total_written -eq 0 ]]; then
    echo "no code symbols matching '$substring' found in $binary" >&2
    exit 1
fi
