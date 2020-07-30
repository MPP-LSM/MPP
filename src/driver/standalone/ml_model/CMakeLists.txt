set(ML_MODEL_PROBLEM_SOURCES
  ml_model_problem.F90
  ml_model_global_vars.F90
  ml_model_meshes.F90
  ml_model_utils.F90
  photosynthesis.F90
  lbl.F90
  lwv.F90
  swv.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/ge)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/soe)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/mpp)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_library)
add_mpp_library(ml_model_problems ${ML_MODEL_PROBLEM_SOURCES})
set_target_properties (ml_model_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES ml_model_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES ml_model_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(ml_model ml_model.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ml_model          DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

