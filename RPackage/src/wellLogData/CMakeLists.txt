cmake_minimum_required(VERSION 3.1)
set(CppFiles register.cpp basicFilterRPackage.cpp fearnheadFilterRPackage.cpp withoutReplacementRPackage.cpp withoutReplacementWithVarianceRPackage.cpp)
set(HeaderFiles basicFilterRPackage.h fearnheadFilterRPackage.h withoutReplacementRPackage.h withoutReplacementWithVarianceRPackage.h)
add_library(wellLogDataRPackage SHARED ${CppFiles} ${HeaderFiles})
target_link_libraries(wellLogDataRPackage PUBLIC wellLogDataCommon Rcpp)
set_target_properties(wellLogDataRPackage PROPERTIES OUTPUT_NAME "wellLogData")
target_compile_definitions(wellLogDataRPackage PRIVATE CUSTOM_STATIC_RCPP)

function(copyBinary buildType destination resultName)
	set(COPIED_BINARY_FILES)
	#Now for MPFR
	get_target_property(MPFR_LOCATION mpfr "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${MPFR_LOCATION}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${MPFR_LOCATION}" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#Now for MPIR
	get_target_property(MPIR_LOCATION mpir "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${MPIR_LOCATION}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${MPIR_LOCATION}" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")
	set(${resultName} ${COPIED_BINARY_FILES} PARENT_SCOPE)
endfunction()

if(WIN32)
	target_compile_definitions(wellLogDataRPackage PUBLIC _CRT_SECURE_NO_WARNINGS)
	target_link_libraries(wellLogDataRPackage PRIVATE "${CMAKE_BINARY_DIR}/RPackage/R.lib")
	set(BUILD_TYPES Release Debug RelWithDebInfo MinSizeRel)
	set(BOOST_BUILD_TYPES RELEASE DEBUG RELEASE RELEASE)
	set(COPIED_BINARY_FILES)
	foreach(val RANGE 3)
		list(GET BUILD_TYPES ${val} BUILD_TYPE)
		list(GET BOOST_BUILD_TYPES ${val} BOOST_BUILD_TYPE)
		string(TOUPPER "${BUILD_TYPE}" UPPER_BUILD_TYPE)
		copyBinary("${BOOST_BUILD_TYPE}" "${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/inst/libs/x64/" additionalBinaryFiles)
		list(APPEND COPIED_BINARY_FILES ${additionalBinaryFiles})
		if(CMAKE_SIZEOF_VOID_P EQUAL 8)
			set_target_properties(wellLogDataRPackage PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${UPPER_BUILD_TYPE} ${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/inst/libs/x64 PDB_OUTPUT_DIRECTORY_${UPPER_BUILD_TYPE} ${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/inst/libs/x64)
		else()
			set_target_properties(wellLogDataRPackage PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${UPPER_BUILD_TYPE} ${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/inst/libs/i386 PDB_OUTPUT_DIRECTORY_${UPPER_BUILD_TYPE} ${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/inst/libs/i386)
		endif()
	endforeach()
	add_custom_target(copyBinary ALL DEPENDS ${COPIED_BINARY_FILES})
else()
	set_target_properties(wellLogDataRPackage PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/RPackage/src)
	set_target_properties(wellLogDataRPackage PROPERTIES PREFIX "")
endif()
