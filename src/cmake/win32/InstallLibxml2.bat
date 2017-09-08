@echo off

set CURRENT_DIR=%~dp0
call %CURRENT_DIR%SetWinBuildEnv.bat

set LIBXML2_SOURCE_DIR=%1/%2
set LIBXML2_SOURCE_DIR=%LIBXML2_SOURCE_DIR:/=\%

set LIBXML2_INSTALL_DIR=%1/%3
set LIBXML2_INSTALL_DIR=%LIBXML2_INSTALL_DIR:/=\%

set LIBXML2_DEBUG=no
if "%4" == "DEBUG" set LIBXML2_DEBUG=yes

pushd %LIBXML2_SOURCE_DIR%\win32

echo -- Configuring libxml2
cscript configure.js prefix=%LIBXML2_INSTALL_DIR% modules=no http=no ftp=no python=no threads=no schematron=no iconv=no debug=%LIBXML2_DEBUG%
if not errorlevel 0 exit /b %errorlevel%

echo -- Patching Makefile and config.h
echo # ================================================================================ > Makefile_Patched
echo # Patched Makefile for libxml2, specific for Windows, Intel C++ compiler and NMake >> Makefile_Patched
echo # >> MakeFile_Patched
echo # Set compiler to Intel C++ compiler and other minor fixes. >> MakeFile_Patched
echo # >> MakeFile_Patched
echo # May 2016, Santosh Prabhakaran ^<sprabhakar@illumina.com^> >> MakeFile_Patched
echo # >> MakeFile_Patched
echo # ================================================================================ >> Makefile_Patched
echo. >> MakeFile_Patched
echo # Set the compiler, archiver and linker to Intel tools. >> Makefile_Patched
echo CPP = icl.exe /EP >> Makefile_Patched
echo CC = icl.exe >> Makefile_Patched
echo LD = xilink.exe >> Makefile_Patched
echo AR = xilib.exe >> Makefile_Patched
echo. >> Makefile_Patched
echo # ================================= End of Patch ================================= >> Makefile_Patched
echo. >> Makefile_Patched

findstr /V "+ OPT:NOWIN98 cl.exe link.exe lib.exe" Makefile >> Makefile_Patched

findstr /V "_snprintf" %LIBXML2_SOURCE_DIR%\config.h > %LIBXML2_SOURCE_DIR%\config.h_Patched
copy /Y %LIBXML2_SOURCE_DIR%\config.h_Patched %LIBXML2_SOURCE_DIR%\config.h

echo -- Building libxml2
nmake /f Makefile_Patched
if not errorlevel 0 exit /b %errorlevel%

echo -- Installing libxml2
nmake /f Makefile_Patched install
if not errorlevel 0 exit /b %errorlevel%

popd
