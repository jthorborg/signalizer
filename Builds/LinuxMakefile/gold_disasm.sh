#!/usr/bin/env bash
#
# Disassemble every code symbol whose (mangled) name contains a given
# substring, one .asm file per symbol, for tracking codegen changes
# through git over time.
#
# Usage: ./gold_disasm.sh <symbol-substring> [path-to-binary]
#
# For each matching symbol, writes:
#   gold/<symbol-substring>/<mangled-name>.asm
# with the c++filt-demangled name as the first line, followed by the
# output of:
#   objdump -d -S -l -M intel -C --no-show-raw-insn --disassemble=<demangled-name> <binary>
#
# Note: nm matches are filtered to code symbols (T/t/W/w/I/i) since
# --disassemble only produces output for those; data symbols (e.g. local
# statics) that happen to share the substring are skipped.
#
# Note: --disassemble=SYMBOL is always an exact match, never a substring or
# prefix. With -C (demangle) enabled, objdump matches --disassemble= against
# the *demangled* name, not the mangled one — passing the mangled name with
# -C silently matches nothing. The demangled name it expects is exactly
# c++filt's output (same libiberty demangler under the hood), so we derive
# it with c++filt and feed it straight back in. This only works reliably
# because the value is kept in a quoted shell variable end-to-end: demangled
# C++ signatures are full of shell metacharacters (<, >, &, (, ), spaces,
# commas) that word-split/redirect/glob if ever unquoted or typed by hand.

set -euo pipefail

usage() {
    echo "Usage: $0 <symbol-substring> [path-to-binary]" >&2
    exit 1
}

[[ $# -ge 1 && $# -le 2 ]] || usage

substring=$1
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
binary=${2:-"$script_dir/build/Signalizer"}
out_dir="$script_dir/gold/$substring"

[[ -f "$binary" ]] || { echo "error: binary not found: $binary" >&2; exit 1; }

mapfile -t symbols < <(
    nm "$binary" | awk '$2 ~ /^[TtWwIi]$/ { print $3 }' | grep -F -- "$substring" | sort -u
)

if [[ ${#symbols[@]} -eq 0 ]]; then
    echo "no code symbols matching '$substring' found in $binary" >&2
    exit 1
fi

mkdir -p "$out_dir"

# Keep filenames well under typical filesystem limits (255 bytes); some
# heavily-templated symbols (e.g. JUCE's X11 symbol loader) mangle to
# well over a thousand characters.
max_name_len=200

for sym in "${symbols[@]}"; do
    if (( ${#sym} > max_name_len )); then
        hash=$(printf '%s' "$sym" | sha1sum | cut -c1-10)
        fname="${sym:0:$((max_name_len - 11))}_$hash.asm"
    else
        fname="$sym.asm"
    fi
    out_file="$out_dir/$fname"

    demangled=$(c++filt "$sym")

    echo "disassembling $sym -> gold/$substring/$fname"
    {
        echo "; $demangled"
        objdump -d -S -l -M intel -C --no-show-raw-insn --disassemble="$demangled" "$binary"
    } > "$out_file"
done

echo "wrote ${#symbols[@]} file(s) to $out_dir"
