set(MLC_PROBLEM_SOURCES
  mlc_problem.F90
  mlc_conditions.F90
  mlc_global_vars.F90
  mlc_meshes.F90
  mlc_parameters.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/ge)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/mlc)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_library)
add_mpp_library(mlc_problems ${MLC_PROBLEM_SOURCES})
set_target_properties (mlc_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mlc_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mlc_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(mlc         mlc.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mlc          DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

