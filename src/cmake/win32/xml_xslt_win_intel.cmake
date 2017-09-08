set (LIBXML_FOLDER_NAME libxml2-${iSAAC_LIBXML2_VERSION})
set (LIBXML_TARBALL ${LIBXML2_REDIST_DIR}/${LIBXML_FOLDER_NAME}.tar.gz)
set (LIBXML_INTEL_DIR libxml2_intel)

# Extract libxml2 from distribution if not already extracted.
if (NOT EXISTS ${BOOTSTRAP_DIR}/${LIBXML_FOLDER_NAME})
	message(STATUS "Libxml2 not found in ${BOOTSTRAP_DIR}.")
	message(STATUS "Libxml2 will be extracted from ${LIBXML_TARBALL}...")
	message("")
	message(STATUS "Please Wait...")

	execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${LIBXML_TARBALL} WORKING_DIRECTORY ${BOOTSTRAP_DIR} RESULT_VARIABLE Result)

	if (Result EQUAL "0")
		message(STATUS "Successfully extracted libxml2 to ${BOOTSTRAP_DIR}.")
	else ()
		message(FATAL_ERROR "Failed to extract libxml2 from ${LIBXML_TARBALL}.")
	endif ()
endif ()

# Build libxml2 using Intel C++ compiler if not already compiled.
if (NOT EXISTS ${BOOTSTRAP_DIR}/${LIBXML_INTEL_DIR})
	message(STATUS "Libxml2 build with Intel compiler not found in ${BOOTSTRAP_DIR}.")
	message(STATUS "Libxml2 will be built from the distribution...")

	execute_process(COMMAND "${CMAKE_SOURCE_DIR}/cmake/win32/InstallLibxml2.bat" "${BOOTSTRAP_DIR}" "${LIBXML_FOLDER_NAME}" "${LIBXML_INTEL_DIR}" "${BUILD_TYPE_UPPER}" WORKING_DIRECTORY ${LIBXML2_REDIST_DIR} RESULT_VARIABLE Result)

	if (Result EQUAL "0")
		message(STATUS "Successfully built libxml2 with Intel compiler to ${BOOTSTRAP_DIR}/${LIBXML_INTEL_DIR}.")
	else ()
		message(FATAL_ERROR "Failed to build libxml2 with Intel compiler.")
	endif ()
else ()
	message(STATUS "Libxml2 build with Intel compiler found at ${BOOTSTRAP_DIR}/${LIBXML_INTEL_DIR}.")
endif ()

include_directories(BEFORE ${BOOTSTRAP_DIR}/${LIBXML_INTEL_DIR}/include)
link_directories(${BOOTSTRAP_DIR}/${LIBXML_INTEL_DIR}/lib)
set (iSAAC_DEP_LIB ${iSAAC_DEP_LIB} libxml2_a.lib)

set (LIBXSLT_FOLDER_NAME libxslt-${iSAAC_LIBXSLT_VERSION})
set (LIBXSLT_TARBALL ${LIBXSLT_REDIST_DIR}/${LIBXSLT_FOLDER_NAME}.tar.gz)
set (LIBXSLT_INTEL_DIR libxslt_intel)

# Extract libxslt from distribution if not already extracted.
if (NOT EXISTS ${BOOTSTRAP_DIR}/${LIBXSLT_FOLDER_NAME})
	message(STATUS "Libxslt not found in ${BOOTSTRAP_DIR}.")
	message(STATUS "Libxslt will be extracted from ${LIBXSLT_TARBALL}...")
	message("")
	message(STATUS "Please Wait...")

	execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${LIBXSLT_TARBALL} WORKING_DIRECTORY ${BOOTSTRAP_DIR} RESULT_VARIABLE Result)

	if (Result EQUAL "0")
		message(STATUS "Successfully extracted libxslt to ${BOOTSTRAP_DIR}.")
	else ()
		message(FATAL_ERROR "Failed to extract libxslt from ${LIBXSLT_TARBALL}.")
	endif ()
endif ()

# Build libxslt using Intel C++ compiler if not already compiled.
if (NOT EXISTS ${BOOTSTRAP_DIR}/${LIBXSLT_INTEL_DIR})
	message(STATUS "Libxslt build with Intel compiler not found in ${BOOTSTRAP_DIR}.")
	message(STATUS "Libxslt will be built from the distribution...")

	execute_process(COMMAND "${CMAKE_SOURCE_DIR}/cmake/win32/InstallLibxslt.bat"
		"${BOOTSTRAP_DIR}" "${LIBXSLT_FOLDER_NAME}" "${LIBXSLT_INTEL_DIR}" "${BOOTSTRAP_DIR}/${LIBXML_INTEL_DIR}/include" "${BOOTSTRAP_DIR}/${LIBXML_INTEL_DIR}/lib" "${BUILD_TYPE_UPPER}"
		WORKING_DIRECTORY ${LIBXSLT_REDIST_DIR} RESULT_VARIABLE Result)

	if (Result EQUAL "0")
		message(STATUS "Successfully built libxslt with Intel compiler to ${BOOTSTRAP_DIR}/${LIBXSLT_INTEL_DIR}.")
	else ()
		message(FATAL_ERROR "Failed to build libxslt with Intel compiler.")
	endif ()
else ()
	message(STATUS "Libxslt build with Intel compiler found at ${BOOTSTRAP_DIR}/${LIBXSLT_INTEL_DIR}.")
endif ()

include_directories(BEFORE ${BOOTSTRAP_DIR}/${LIBXSLT_INTEL_DIR}/include)
link_directories(${BOOTSTRAP_DIR}/${LIBXSLT_INTEL_DIR}/lib)
set (iSAAC_DEP_LIB ${iSAAC_DEP_LIB} libxslt_a.lib libexslt_a.lib)
