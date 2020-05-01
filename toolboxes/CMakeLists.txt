
if (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
endif (WIN32)

add_subdirectory(core)
add_subdirectory(cmr)