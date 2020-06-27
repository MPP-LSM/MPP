set(MPP_GE_SOURCES
  CanopyTurbulence.F90
  GoveqnCanopyAirTemperatureType.F90
  GoveqnCanopyAirVaporType.F90
  GoveqnCanopyLeafTemperatureType.F90
  GoveqnLeafBoundaryLayer.F90
  GoveqnRichardsODEPressureType.F90
  GoveqnThermalEnthalpySoilType.F90
  GoveqnThermalKSPTemperatureSSWType.F90
  GoveqnThermalKSPTemperatureSnowType.F90
  GoveqnThermalKSPTemperatureSoilType.F90
  MathTools.F90
  RichardsMod.F90
  ThermalEnthalpyMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/auxvar)

include(add_mpp_library)
add_mpp_library(mpp_ge ${MPP_GE_SOURCES})
set_target_properties (mpp_ge PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_ge; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_ge; ${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_ge DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/ge)
endif()

