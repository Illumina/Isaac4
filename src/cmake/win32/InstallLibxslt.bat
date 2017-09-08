@echo off

set CURRENT_DIR=%~dp0
call %CURRENT_DIR%SetWinBuildEnv.bat

set LIBXSLT_SOURCE_DIR=%1/%2
set LIBXSLT_SOURCE_DIR=%LIBXSLT_SOURCE_DIR:/=\%

set LIBXSLT_INSTALL_DIR=%1/%3
set LIBXSLT_INSTALL_DIR=%LIBXSLT_INSTALL_DIR:/=\%

set LIBXML2_INCLUDE_DIR=%4
set LIBXML2_INCLUDE_DIR=%LIBXML2_INCLUDE_DIR:/=\%

set LIBXML2_LIB_DIR=%5
set LIBXML2_LIB_DIR=%LIBXML2_LIB_DIR:/=\%

set LIBXSLT_DEBUG=no
if "%6" == "DEBUG" set LIBXSLT_DEBUG=yes

pushd %LIBXSLT_SOURCE_DIR%\win32

echo -- Configuring libxslt
cscript configure.js prefix=%LIBXSLT_INSTALL_DIR% include=%LIBXML2_INCLUDE_DIR% lib=%LIBXML2_LIB_DIR% crypto=no iconv=no debug=%LIBXSLT_DEBUG%
if not errorlevel 0 exit /b %errorlevel%

echo -- Patching Makefile and config.h
echo # ================================================================================ > Makefile_Patched
echo # Patched Makefile for libxslt, specific for Windows, Intel C++ compiler and NMake >> Makefile_Patched
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

findstr /V "OPT:NOWIN98 cl.exe link.exe lib.exe" Makefile >> Makefile_Patched

findstr /V "_snprintf" %LIBXSLT_SOURCE_DIR%\libxslt\win32config.h > %LIBXSLT_SOURCE_DIR%\libxslt\win32config.h_Patched
copy /Y %LIBXSLT_SOURCE_DIR%\libxslt\win32config.h_Patched %LIBXSLT_SOURCE_DIR%\libxslt\win32config.h

echo -- Building libxslt
nmake /f Makefile_Patched
if not errorlevel 0 exit /b %errorlevel%

echo -- Installing libxslt
nmake /f Makefile_Patched install
if not errorlevel 0 exit /b %errorlevel%

popd
