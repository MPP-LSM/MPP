set(MPP_SOE_SOURCES
  SystemOfEquationsBasePointerType.F90
  SystemOfEquationsBaseType.F90
  SystemOfEquationsMLCAuxType.F90
  SystemOfEquationsMLCType.F90
  SystemOfEquationsLblAuxType.F90
  SystemOfEquationsLblType.F90
  SystemOfEquationsPhotosynthesisAuxType.F90
  SystemOfEquationsPhotosynthesisType.F90
  SystemOfEquationsThermalEnthalpyType.F90
  SystemOfEquationsThermalAuxMod.F90
  SystemOfEquationsThermalType.F90
  SystemOfEquationsTHType.F90
  SystemOfEquationsVSFMType.F90  
  SystemOfEquationsThermalEnthalpyAuxMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/auxvar)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/ge)

include(add_mpp_library)
add_mpp_library(mpp_soe ${MPP_SOE_SOURCES})
set_target_properties (mpp_soe PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_soe; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_soe; ${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_soe DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/soe)
endif()

