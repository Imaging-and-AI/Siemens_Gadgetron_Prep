# this toolbox contains some key algorithms for cmr related applications

if ( HAS_64_BIT )

    if (WIN32)
        ADD_DEFINITIONS(-D__BUILD_GADGETRON_SIEMENS_TOOLBOX_CMR__)
    endif (WIN32)

    if(WIN32)
        link_directories(${Boost_LIBRARY_DIRS})
    endif(WIN32)

    include_directories(
        ${ACE_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        ${ISMRMRD_INCLUDE_DIR}
        ${FFTW3_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
        ${MKL_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/toolboxes
        ${CMAKE_SOURCE_DIR}/toolboxes/cmr
        ${CMAKE_SOURCE_DIR}/gadgets
    )

    INCLUDE_DIRECTORIES(${GADGETRON_INCLUDE_DIR})

    set( cmr_header_files gadgetron_siemens_toolbox_cmr_export.h )
    set( cmr_src_files )

    set( application_header_files
        application/cmr_sasha_t1_t2_mapping.h
        application/cmr_multi_parametric_mapping.h
        )

    set( application_src_files 
        application/cmr_sasha_t1_t2_mapping.cpp
        application/cmr_multi_parametric_mapping.cpp
        )

    set( schema_files )

    source_group(schema FILES ${schema_files})

    set(source_files ${cmr_header_files}
            ${cmr_src_files}
            ${application_header_files}
            ${application_src_files} )

    source_group(application FILES  ${application_header_files} ${application_src_files})

    add_library(gadgetron_siemens_toolbox_cmr SHARED ${source_files})

    set_target_properties(gadgetron_siemens_toolbox_cmr PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

    target_link_libraries(gadgetron_siemens_toolbox_cmr
                          ${MKL_LIBRARIES}
                          ${MATLAB_LIBRARIES}
                          ${Boost_LIBRARIES}
                          ${ISMRMRD_LIBRARIES}
                          optimized  gadgetron_toolbox_cpufft debug gadgetron_toolbox_cpufft${CMAKE_DEBUG_SUFFIX}
                          optimized  gadgetron_toolbox_cpuklt debug gadgetron_toolbox_cpuklt${CMAKE_DEBUG_SUFFIX}
                          optimized  gadgetron_toolbox_image_analyze_io debug gadgetron_toolbox_image_analyze_io${CMAKE_DEBUG_SUFFIX}
                          optimized gadgetron_toolbox_cpucore debug gadgetron_toolbox_cpucore${CMAKE_DEBUG_SUFFIX}
                          optimized gadgetron_toolbox_cpucore_math debug gadgetron_toolbox_cpucore_math${CMAKE_DEBUG_SUFFIX}
                          optimized gadgetron_toolbox_mri_core debug gadgetron_toolbox_mri_core${CMAKE_DEBUG_SUFFIX}
                          optimized gadgetron_toolbox_log debug gadgetron_toolbox_log${CMAKE_DEBUG_SUFFIX}
                          optimized gadgetron_toolbox_cpudwt debug gadgetron_toolbox_cpudwt${CMAKE_DEBUG_SUFFIX} 
                          optimized gadgetron_toolbox_cmr debug gadgetron_toolbox_cmr${CMAKE_DEBUG_SUFFIX} 
                          optimized gadgetron_toolbox_pr debug gadgetron_toolbox_pr${CMAKE_DEBUG_SUFFIX} 
                          optimized gadgetron_toolbox_cpu_image debug gadgetron_toolbox_cpu_image${CMAKE_DEBUG_SUFFIX} )

    install(TARGETS gadgetron_siemens_toolbox_cmr DESTINATION lib)

endif ( HAS_64_BIT )
