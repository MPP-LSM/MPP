set(VSFM_DRIVER_SOURCES
  vsfm_celia1990.F90
  vsfm_sy1991.F90
  vsfm_vchannel.F90
  vsfm_vchannel_operator_split.F90
  vsfm_manoli2014.F90
  )

set(VSFM_PROBLEM_SOURCES
  vsfm_celia1990_problem.F90
  vsfm_sy1991_problem.F90
  vsfm_vchannel_problem.F90
  vsfm_vchannel_problem_operator_split.F90
  vsfm_manoli2014_problem.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_library)
add_mpp_library(vsfm_problems ${VSFM_PROBLEM_SOURCES})
set_target_properties (vsfm_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES vsfm_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES vsfm_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(vsfm_celia1990 vsfm_celia1990.F90)
add_mpp_executable(vsfm_sy1991 vsfm_sy1991.F90)
add_mpp_executable(vsfm_vchannel vsfm_vchannel.F90)
add_mpp_executable(vsfm_vchannel_operator_split vsfm_vchannel_operator_split.F90)
add_mpp_executable(vsfm_manoli2014 vsfm_manoli2014.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS vsfm_problems DESTINATION lib)

  install(TARGETS vsfm_celia1990 DESTINATION bin)
  install(TARGETS vsfm_sy1991 DESTINATION bin)
  install(TARGETS vsfm_vchannel DESTINATION bin)
  install(TARGETS vsfm_vchannel_operator_split DESTINATION bin)
  install(TARGETS vsfm_manoli2014 DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

