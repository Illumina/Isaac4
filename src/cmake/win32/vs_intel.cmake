#
# Compiler definitions and flags for Visual Studio with Intel C++ Compiler.
#
add_definitions (-DHAVE_ERFCF)
add_compile_options (/EHsc /Wall)
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:8388608")

if (BUILD_TYPE_UPPER MATCHES DEBUG)
	add_compile_options (/ZI /Od /Ob2)
endif ()

configure_file(${CMAKE_SOURCE_DIR}/cmake/win32/Make.bat.tmpl ${CMAKE_CURRENT_BINARY_DIR}/Make.bat)
