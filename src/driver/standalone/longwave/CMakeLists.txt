set(LONGWAVE_PROBLEM_SOURCES
  longwave_problem.F90
  longwave_global_vars.F90
  longwave_meshes.F90
  longwave_parameters.F90
  longwave_conditions.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/ge)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/soe)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/mpp)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/util)

link_libraries(${MPP_LIBRARIES})
link_directories(${PETSC_LIBRARY_DIRS})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_library)
add_mpp_library(LONGWAVE_PROBLEMs ${LONGWAVE_PROBLEM_SOURCES})
set_target_properties (LONGWAVE_PROBLEMs PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES LONGWAVE_PROBLEMs; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES LONGWAVE_PROBLEMs; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(longwave longwave.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS longwave DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

