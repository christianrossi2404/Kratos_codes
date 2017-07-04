set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosHeatEquationApplication ****")
message("HEAT_EQUATION_APPLICATION.....................${HEAT_EQUATION_APPLICATION}")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_HEATEQUATION_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/heat_equation_application.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/heat_equation_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/heat_equation.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element_2d.h

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/heat_equation_Neumann_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosHeatEquationApplication SHARED ${KRATOS_HEATEQUATION_APPLICATION_SOURCES})
target_link_libraries(KratosHeatEquationApplication KratosCore )
set_target_properties(KratosHeatEquationApplication PROPERTIES PREFIX "")
install(TARGETS KratosHeatEquationApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosHeatEquationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

#if(${HEAT_EQUATION_APPLICATION} MATCHES ON)
#   add_subdirectory(heat_equation_application)
#endif(${HEAT_EQUATION_APPLICATION} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/HeatEquationApplication.py" DESTINATION KratosMultiphysics )
