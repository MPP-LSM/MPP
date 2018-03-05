set(THERMAL_PROBLEM_SOURCES
  thermal_mms_problem.F90
  thermal_mms_vars.F90
  thermal_mms_steady_state_problem_1D.F90
  thermal_mms_steady_state_problem_2D.F90
  thermal_mms_steady_state_problem_3D.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/th)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal-e)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_library)
add_mpp_library(thermal_problems ${THERMAL_PROBLEM_SOURCES})
set_target_properties (thermal_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES thermal_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES thermal_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(thermal_mms thermal_mms.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS thermal_problems DESTINATION lib)

  install(TARGETS thermal_mms DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

