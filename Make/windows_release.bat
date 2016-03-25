@echo off

echo Building Signalizer v. X release targets

REM set up x86 target
call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86
devenv "..\builds\visualstudio2010\signalizer.sln" /build "Release|win32" /LINK:/VERSION:major.minor

REM set up x64 target
call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x64
devenv "..\builds\visualstudio2010\signalizer.sln" /build "Release|x64" /LINK:/VERSION:major.minor

echo.
echo All builds completed