set(LEAFBNDLYR_PROBLEM_SOURCES
  leafbndlyr_problem.F90
  lbl_global_vars.F90
  lbl_meshes.F90
  lbl_parameters.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/ge)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/soe)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/mpp)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_library)
add_mpp_library(leafbndlyr_problems ${LEAFBNDLYR_PROBLEM_SOURCES})
set_target_properties (leafbndlyr_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES leafbndlyr_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES leafbndlyr_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(leafbndlyr leafbndlyr.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS leafbndlyr DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

