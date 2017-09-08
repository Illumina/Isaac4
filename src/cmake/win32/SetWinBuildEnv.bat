@echo off

REM Only this file needs to be modified if any version change in the Visual Studio or Intel compiler.
REM Also make sure to delete the contents of the %TEMP% folder because Boost.Build caches files here.

REM If this batch file was already invoked, simply exit.
if not "%BOOST_TOOLSET%" == "" exit /b 0

REM Set Visual Studio environment.
call "%ProgramFiles(x86)%\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" AMD64

REM Set Intel compiler environment.
call "%ICPP_COMPILER16%bin\ipsxe-comp-vars.bat" intel64 vs2015

set BOOST_TOOLSET=intel-16.0-vc14
