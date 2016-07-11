add_subdirectory(util)
add_subdirectory(dtypes)
add_subdirectory(vsfm)
add_subdirectory(thermal)
add_subdirectory(thermal-e)

set(MPP_LIBRARIES mpp_dtypes; mpp_thermal; mpp_thermal_e; mpp_util; mpp_vsfm; ${MPP_LIBRARIES} PARENT_SCOPE)

