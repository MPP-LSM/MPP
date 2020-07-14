set(MPP_AUXVAR_SOURCES
  CanopyAirTemperatureAuxType.F90
  CanopyAirTemperatureConnAuxType.F90
  CanopyAirVaporAuxType.F90
  CanopyAirVaporConnAuxType.F90
  CanopyLeafTemperatureAuxType.F90
  CanopyTurbulenceAuxType.F90
  LeafBoundaryLayerAuxType.F90
  PhotosynthesisAuxVarType.F90
  RichardsODEPressureAuxType.F90
  RichardsODEPressureAuxMod.F90
  RichardsODEPressureConnAuxType.F90
  SystemOfEquationsVSFMAuxType.F90
  SystemOfEquationsThermalEnthalpyAuxType.F90
  SystemOfEquationsThermalAuxType.F90
  ThermalEnthalpySoilAuxType.F90
  ThermalEnthalpySoilAuxMod.F90
  ThermalKSPTemperatureSSWAuxMod.F90
  ThermalKSPTemperatureBaseAuxType.F90
  ThermalKSPTemperatureSSWAuxType.F90
  ThermalKSPTemperatureSnowAuxMod.F90
  ThermalKSPTemperatureSnowAuxType.F90
  ThermalKSPTemperatureSoilAuxMod.F90
  ThermalKSPTemperatureSoilAuxType.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)

include(add_mpp_library)
add_mpp_library(mpp_auxvar ${MPP_AUXVAR_SOURCES})
set_target_properties (mpp_auxvar PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_auxvar;${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_auxvar;${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_auxvar DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/vsfm)
endif()

