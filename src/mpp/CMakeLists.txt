add_subdirectory(util)
add_subdirectory(dtypes)
add_subdirectory(auxvar)
add_subdirectory(vsfm)
add_subdirectory(thermal)
add_subdirectory(thermal-e)
add_subdirectory(th)
add_subdirectory(mlc)

set(MPP_LIBRARIES mpp_dtypes; mpp_auxvar; mpp_thermal; mpp_thermal_e; mpp_util; mpp_vsfm; mpp_th; mpp_mlc; ${MPP_LIBRARIES} PARENT_SCOPE)

