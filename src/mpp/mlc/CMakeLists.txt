set(MPP_MLC_SOURCES
  CanopyTurbulence.F90
  GoveqnCanopyAirTemperatureType.F90
  GoveqnCanopyAirVaporType.F90
  GoveqnCanopyLeafTemperatureType.F90
  MultiPhysicsProbMLC.F90
  SystemOfEquationsMLCType.F90
  MathTools.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/auxvar)

include(add_mpp_library)
add_mpp_library(mpp_mlc ${MPP_MLC_SOURCES})
set_target_properties (mpp_mlc PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_mlc; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_mlc; ${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_mlc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/mlc)
endif()

