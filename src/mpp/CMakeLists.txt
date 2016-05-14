add_subdirectory(util)
add_subdirectory(dtypes)
add_subdirectory(thermal)
add_subdirectory(vsfm)

set(MPP_LIBRARIES mpp_dtypes; mpp_thermal; mpp_util; mpp_vsfm; ${MPP_LIBRARIES} PARENT_SCOPE)

