set(PHOTOSYNTHESIS_PROBLEM_SOURCES
  photosynthesis_problem.F90
  photosynthesis_global_vars.F90
  photosynthesis_meshes.F90
  photosynthesis_parameters.F90
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
add_mpp_library(photosynthesis_problems ${PHOTOSYNTHESIS_PROBLEM_SOURCES})
set_target_properties (photosynthesis_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES photosynthesis_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES photosynthesis_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(photosynthesis photosynthesis.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS photosynthesis DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

