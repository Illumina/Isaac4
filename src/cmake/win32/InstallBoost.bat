@echo off

set CURRENT_DIR=%~dp0
call %CURRENT_DIR%SetWinBuildEnv.bat

set BOOST_SOURCE_DIR=%1/%2
set BOOST_SOURCE_DIR=%BOOST_SOURCE_DIR:/=\%

set BOOST_INSTALL_DIR=%1/%3
set BOOST_INSTALL_DIR=%BOOST_INSTALL_DIR:/=\%

set ZLIB_INSTALL_DIR=%1/%4
set ZLIB_INSTALL_DIR=%ZLIB_INSTALL_DIR:/=\%

set BOOST_VARIANT=release
if "%5" == "DEBUG" set BOOST_VARIANT=debug

pushd %BOOST_SOURCE_DIR%

call bootstrap.bat

set iSAAC_BOOST_LIBRARIES=
for %%L in (%iSAAC_BOOST_COMPONENTS%) do call :concat --with-%%L

b2 --build-dir=%BOOST_INSTALL_DIR% --prefix=%BOOST_INSTALL_DIR% toolset=%BOOST_TOOLSET% variant=%BOOST_VARIANT% threading=multi architecture=x86 address-model=64 define=_WIN32_WINNT=0x0601 install %iSAAC_BOOST_LIBRARIES% -s ZLIB_SOURCE=%ZLIB_INSTALL_DIR%

popd
goto :eof

REM Create command line arguments for b2 in the format [--with-<library> --with-<library> ...] where the list of libraries are taken from the iSAAC_BOOST_COMPONENTS variable.
:concat
set iSAAC_BOOST_LIBRARIES=%iSAAC_BOOST_LIBRARIES% %1
goto :eof
