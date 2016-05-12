add_subdirectory(util)
add_subdirectory(dtypes)
add_subdirectory(thermal)
add_subdirectory(vsfm)

set(MPP_LIBRARIES dtypes; thermal; util; vsfm; ${MPP_LIBRARIES} PARENT_SCOPE)

