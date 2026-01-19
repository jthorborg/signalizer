# Signalizer

Public repository for the real-time audio visualization plugin Signalizer.

More info can be found at this page: www.jthorborg.com/index.html?ipage=signalizer

Pre-built binaries can be found at this page: https://github.com/jthorborg/signalizer/releases

![signalizer_screensaver.png](https://bitbucket.org/repo/jnBRk8/images/675350869-signalizer_screensaver.png)

## Building Signalizer

The currently supported build platforms are Windows and OS X. If you haven't yet, run:
`$ python3 prepare.py`

Then you can run:
`$ python3 build.py -[hdsjnb]`

And a zipped Signalizer release for your platform will be built into `/Releases`.
You can also use either of the platform specific solutions in `/Builds/` to do development.
The `*.jucer` file is used to rebuild the solutions, and the python build system invokes the solutions and packages releases.

If you wish to compile VST 2s, you will need to acquire the SDK and place it in `../SDKs/vstsdk2.4`. 
Similarly for other proprietary platforms.

Linux support is experimental and only confirmed to work on Ubuntu 24.
On macOS, if you want to use the build scripts you need to have the Xcode command line tools installed.
For Windows, you will need Visual Studio 2022+.

## Performance

## Known issues

- "Plugin is damaged and can't be opened" on newer versions of macOS (due to lack of notarization): See the commentary in [this file](Make/macos_installation_advice.txt)