call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat" %1
devenv "..\builds\visualstudio2010\signalizer.sln" /build %2

exit %ERRORLEVEL%
