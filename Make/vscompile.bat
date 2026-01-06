call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" %1
devenv "..\Builds\VisualStudio2022\Signalizer.sln" /build %2

exit %ERRORLEVEL%
