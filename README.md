# Signalizer

Public repository for the real-time audio visualization plugin Signalizer.

More info can be found at this page: www.jthorborg.com/index.html?ipage=signalizer

Pre-built binaries can be found at this page: https://github.com/jthorborg/signalizer/releases

![signalizer_screensaver.png](https://bitbucket.org/repo/jnBRk8/images/675350869-signalizer_screensaver.png)

## Building Signalizer

First-time setup:

```
python prepare.py
```

**Development build** (Debug Standalone, fastest - good for iteration and testing):

```
python build.py dev
```

**Release build** (all plugin formats, packaged into `/Releases`):

```
python build.py release
```

Both commands accept `--verbose` to stream the full compiler output. Run `python3 build.py dev --help` or `release --help` for all options.

Platform notes:
- Windows: requires Visual Studio 2022
- macOS: requires Xcode command line tools
- Linux: experimental, confirmed on Ubuntu 24

VST2 requires the SDK at `../SDKs/vstsdk2.4`.

See more repository / workflow information in [AGENTS.md](AGENTS.md).

## Performance

## Known issues

- "Plugin is damaged and can't be opened" on newer versions of macOS (due to lack of notarization): See the commentary in [this file](Make/macos_installation_advice.txt)