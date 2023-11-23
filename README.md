# Signalizer

Public repository for the real-time audio visualization plugin Signalizer.
More info can be found at this page: www.jthorborg.com/index.html?ipage=signalizer
Pre-built binaries can be found at this page: https://bitbucket.org/Mayae/signalizer/downloads/

![signalizer_screensaver.png](https://bitbucket.org/repo/jnBRk8/images/675350869-signalizer_screensaver.png)

## Building Signalizer

The currently supported build platforms are Windows and OS X. If you haven't yet, run:
`$ python prepare.py`

Then from `/Make`, run:
`$ python build_[win|linux|osx].py [-inc:major|minor|patch]`

And a zipped Signalizer release will be built into `/Releases`.

Linux support is experimental and only confirmed to work on Ubuntu 16. See Make/LinuxInstructions.txt for instructions.
For OS X, you will need Xcode 5.1.1+ as well as the audio unit SDKs installed.
If you want to use the build scripts you need to have the Xcode command line tools installed.

If you wish to compile VSTs, you will need to acquire the VST3 SDK. Similarly for other proprietary platforms.
For Windows, you will need Visual Studio 2022+.

Additionally, you may need to set up include directories to point to your specific SDK locations.

## Performance

Program maxing out a whole core with VSync enabled? If you are using NVidia graphics, it is a driver issue.. Or, "feature", as they call it. Check this out: http://www.retrocopy.com/blog/29/nvidia-threaded-optimization-oxymoron.aspx
Notice, that it is not actually using all of the CPU, it is merely busywaiting (but yielding). What this means in reality is, that your core will be maxed out, but everything will remain as responsive. See this: http://forum.openscenegraph.org/viewtopic.php?t=3653#18283

## Known issues
