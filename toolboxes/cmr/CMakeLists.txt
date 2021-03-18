# this toolbox contains some key algorithms for cmr related applications

include_directories(
    ${GADGETRON_INCLUDE_DIR}
    )

include_directories(${CMAKE_SOURCE_DIR}
                ${ISMRMRD_INCLUDE_DIR}
                ${Boost_INCLUDE_DIR}
                ${PYTHON_INCLUDE_PATH}
                ${NUMPY_INCLUDE_DIRS}
                ${ARMADILLO_LIBRARIES} 
                )

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
                        gadgetron_core
                        gadgetron_toolbox_cpufft
                        gadgetron_toolbox_cpuklt
                        gadgetron_toolbox_image_analyze_io
                        gadgetron_toolbox_cpucore
                        gadgetron_toolbox_cpucore_math
                        gadgetron_toolbox_mri_core
                        gadgetron_toolbox_log
                        gadgetron_toolbox_cpudwt
                        gadgetron_toolbox_cmr
                        gadgetron_toolbox_pr
                        gadgetron_toolbox_cpu_image 
                        ${ISMRMRD_LIBRARIES}
                        Boost::system
                        Boost::filesystem
                        Boost::program_options
                        )

target_include_directories(gadgetron_siemens_toolbox_cmr PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/application>
        )

install(TARGETS gadgetron_siemens_toolbox_cmr DESTINATION lib)
