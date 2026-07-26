# Signalizer

Public repository for the real-time audio visualization plugin Signalizer.
Download prebuilt releases of Signalizer here: https://github.com/jthorborg/signalizer/releases

More info can be found at this page: www.jthorborg.com/index.html?ipage=signalizer

![signalizer_screensaver.png](https://bitbucket.org/repo/jnBRk8/images/675350869-signalizer_screensaver.png)

## Building Signalizer

First-time setup:

```
python prepare.py
```

**Development build** (Debug Standalone, fastest - good for iteration and testing):

```
python build.py dev [--optimized] [--run]
```

**Release build** (all plugin formats, packaged into `/Releases`):

```
python build.py release
```

Both commands accept `--verbose` to stream the full compiler output. Run `python3 build.py dev --help` or `release --help` for all options.

Platform notes:
- Windows: requires Visual Studio 2022, targets Windows 10+
- macOS: requires Xcode command line tools, targets 10.13+
- Linux: any apt-based distribution. Requires GCC 12 or newer; `prepare.py` installs it when the distribution default is older. Verified on Ubuntu 22.04 and 24.04.

VST2 requires the SDK at `../SDKs/vstsdk2.4`.

### Linux binary compatibility

The build host bounds the compatibility of the resulting binaries.

| Built on | Requires | Covers |
| --- | --- | --- |
| Ubuntu 22.04 | glibc 2.35, `GLIBCXX_3.4.30` | Ubuntu 22.04+, Debian 12+, Fedora 36+, and derivatives |
| Ubuntu 24.04 | glibc 2.38, `GLIBCXX_3.4.32` | Ubuntu 23.10+, Debian 13+, Fedora 39+ |

Verify what a given build actually requires with `readelf -V <binary>`.

See more repository / workflow information in [AGENTS.md](AGENTS.md).

## Performance

## Known issues

- "Plugin is damaged and can't be opened" on newer versions of macOS (due to lack of notarization): See the commentary in [this file](Make/macos_installation_advice.txt)
