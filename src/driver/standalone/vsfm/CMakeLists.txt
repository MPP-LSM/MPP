set(VSFM_DRIVER_SOURCES
  vsfm_celia1990.F90
  vsfm_vchannel.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)

#include(add_mpp_library)
#add_mpp_library(vsfm_driver ${VSFM_DRIVER_SOURCES})

#set(MPP_LIBRARIES vsfm_driver;${MPP_LIBRARIES} PARENT_SCOPE)
#set(MPP_LIBRARIES vsfm_driver;${MPP_LIBRARIES})

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(vsfm_celia1990 vsfm_celia1990.F90)
add_mpp_executable(vsfm_vchannel vsfm_vchannel.F90)

message ("MPP_LIBRARIES = ${MPP_LIBRARIES} ")

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS vsfm_celia1990 DESTINATION bin)
  install(TARGETS vsfm_vchannel DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

