if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETS__)
  link_directories(${Boost_LIBRARY_DIRS})
endif (WIN32)

# These dependencies should have been found already to even get in here, but we play safe
find_package(ACE)
find_package(FFTW3 COMPONENTS single double)

include_directories(
    ${ACE_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
  )

add_subdirectory(sasha_hc)
add_subdirectory(moco_ave)