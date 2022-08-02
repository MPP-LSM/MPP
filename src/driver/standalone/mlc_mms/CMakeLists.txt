set(MLC_MMS_PROBLEM_SOURCES
  mlc_mms_problem.F90
  mlc_mms_conditions.F90
  mlc_mms_global_vars.F90
  mlc_mms_meshes.F90
  mlc_mms_parameters.F90
  mlc_mms_functions.F90
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
add_mpp_library(mlc_mms_problems ${MLC_MMS_PROBLEM_SOURCES})
set_target_properties (mlc_mms_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mlc_mms_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mlc_mms_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(mlc_mms         mlc_mms.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mlc_mms          DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

