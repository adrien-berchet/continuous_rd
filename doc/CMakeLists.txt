FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)

	# Copy the Doxygen config file and replaces the variables of the form @VAR@ but not the variables of the form ${VAR}
	CONFIGURE_FILE(
		${CMAKE_CURRENT_SOURCE_DIR}/doxygen.conf.in
		${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
		@ONLY
	)

	# Copy the README file that Doxygen uses as mainpage.
	CONFIGURE_FILE(
		${CMAKE_CURRENT_SOURCE_DIR}/../README.md
		${CMAKE_CURRENT_BINARY_DIR}/README.md
		COPYONLY
	)

	# Custom command used to build Doxygen documentation.
	ADD_CUSTOM_COMMAND(
		OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
		COMMAND ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
		COMMENT "Generating Doxygen html output..."
		DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
	)

	# New target in Makefile that can be executed with "make doc"
	ADD_CUSTOM_TARGET( doc ALL
		COMMAND "echo"
		COMMENT "Documentation was generated here: ${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
		DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
	)

	INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ DESTINATION doc )

endif(DOXYGEN_FOUND)
