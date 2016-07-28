set(VSFM_DRIVER_SOURCES
  heat_transport_1D.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/th)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal-e)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(heat_transport_1D heat_transport_1D.F90)
add_mpp_executable(heat_transport_1D_with_advection heat_transport_1D_with_advection.F90)
add_mpp_executable(mass_and_heat_model mass_and_heat_model.F90)

message ("MPP_LIBRARIES = ${MPP_LIBRARIES} ")

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS heat_transport_1D DESTINATION bin)
  install(TARGETS heat_transport_1D_with_advection DESTINATION bin)
  install(TARGETS mass_and_heat_model DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()
