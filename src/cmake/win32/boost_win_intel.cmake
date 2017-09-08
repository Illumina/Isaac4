set (BOOTSTRAP_DIR ${CMAKE_CURRENT_BINARY_DIR}/bootstrap)
file (MAKE_DIRECTORY ${BOOTSTRAP_DIR})

set (ZLIB_FOLDER_NAME zlib-${iSAAC_ZLIB_VERSION})
set (ZLIB_TARBALL ${ZLIB_REDIST_DIR}/${ZLIB_FOLDER_NAME}.tar.gz)

# Extract zlib from distribution if not already extracted.
if (NOT EXISTS ${BOOTSTRAP_DIR}/${ZLIB_FOLDER_NAME})
	message(STATUS "ZLib not found in ${BOOTSTRAP_DIR}.")
	message(STATUS "ZLib will be extracted from ${ZLIB_TARBALL}...")
	message("")
	message(STATUS "Please Wait...")

	execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${ZLIB_TARBALL} WORKING_DIRECTORY ${BOOTSTRAP_DIR} RESULT_VARIABLE Result)

	if (Result EQUAL "0")
		message(STATUS "Successfully extracted zlib to ${BOOTSTRAP_DIR}.")
	else ()
		message(FATAL_ERROR "Failed to extract zlib from ${ZLIB_TARBALL}.")
	endif ()
endif ()

include_directories(BEFORE ${BOOTSTRAP_DIR}/${ZLIB_FOLDER_NAME})

string (REGEX REPLACE "\\." "_" BOOST_VERSION ${iSAAC_BOOST_VERSION})

set (BOOST_FOLDER_NAME boost_${BOOST_VERSION})
set (BOOST_TARBALL ${BOOST_REDIST_DIR}/${BOOST_FOLDER_NAME}.tar.bz2)
set (BOOST_INTEL_DIR boost_intel)

# Extract boost from distribution if not already extracted.
if (NOT EXISTS ${BOOTSTRAP_DIR}/${BOOST_FOLDER_NAME})
	message(STATUS "Boost not found in ${BOOTSTRAP_DIR}.")
	message(STATUS "Boost will be extracted from ${BOOST_TARBALL}...")
	message("")
	message(STATUS "Please Wait...")

	execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${BOOST_TARBALL} WORKING_DIRECTORY ${BOOTSTRAP_DIR} RESULT_VARIABLE Result)

	if (Result EQUAL "0")
		message(STATUS "Successfully extracted boost to ${BOOTSTRAP_DIR}.")
	else ()
		message(FATAL_ERROR "Failed to extract boost from ${BOOST_TARBALL}.")
	endif ()
endif ()

# Build boost using Intel C++ compiler if not already compiled.
if (NOT EXISTS ${BOOTSTRAP_DIR}/${BOOST_INTEL_DIR})
	message(STATUS "Boost build with Intel compiler not found in ${BOOTSTRAP_DIR}.")
	message(STATUS "Boost will be built from the distribution...")

    set(ENV{iSAAC_BOOST_COMPONENTS} "${iSAAC_BOOST_COMPONENTS}")
	execute_process(COMMAND "${CMAKE_SOURCE_DIR}/cmake/win32/InstallBoost.bat" "${BOOTSTRAP_DIR}" "${BOOST_FOLDER_NAME}" "${BOOST_INTEL_DIR}" "${ZLIB_FOLDER_NAME}" "${BUILD_TYPE_UPPER}" WORKING_DIRECTORY ${BOOST_REDIST_DIR} RESULT_VARIABLE Result)

	if (Result EQUAL "0")
		message(STATUS "Successfully built boost with Intel compiler to ${BOOTSTRAP_DIR}/${BOOST_INTEL_DIR}.")
	    set (BOOST_ROOT "${BOOTSTRAP_DIR}/${BOOST_INTEL_DIR}")
	else ()
		message(FATAL_ERROR "Failed to build boost with Intel compiler.")
	endif ()
else ()
	message(STATUS "Boost build with Intel compiler found at ${BOOTSTRAP_DIR}/${BOOST_INTEL_DIR}.")
endif ()

string (FIND "${BOOST_VERSION}" "_" US_POS REVERSE)
string (SUBSTRING "${BOOST_VERSION}" 0 ${US_POS}-1 BOOST_MAJOR_VERSION)

set (BOOST_VERSION_DIR boost-${BOOST_MAJOR_VERSION})

# Patch vector.hpp in the boost headers (Compiler is not able to resolve value_at_c type to the correct namespace)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${BOOST_REDIST_DIR}/patch/vector.hpp ${BOOTSTRAP_DIR}/${BOOST_INTEL_DIR}/include/${BOOST_VERSION_DIR}/boost/fusion/container/vector)

include_directories(BEFORE ${BOOTSTRAP_DIR}/${BOOST_INTEL_DIR}/include/${BOOST_VERSION_DIR})
link_directories(${BOOTSTRAP_DIR}/${BOOST_INTEL_DIR}/lib)
