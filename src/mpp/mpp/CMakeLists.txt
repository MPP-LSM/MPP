set(MPP_MPP_SOURCES
  MultiPhysicsProbMLC.F90
  MultiPhysicsProbTH.F90
  MultiPhysicsProbThermal.F90
  MultiPhysicsProbThermalEnthalpy.F90
  MultiPhysicsProbVSFM.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/auxvar)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/ge)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/soe)

include(add_mpp_library)
add_mpp_library(mpp_mpp ${MPP_MPP_SOURCES})
set_target_properties (mpp_mpp PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_mpp; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_mpp; ${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_mpp DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/mpp)
endif()

