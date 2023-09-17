call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" %1
devenv "..\builds\visualstudio2010\signalizer.sln" /build %2

exit %ERRORLEVEL%
