# this toolbox contains some key algorithms for cmr related applications

if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_SIEMENS_TOOLBOX_CORE__)
endif (WIN32)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
    ${GADGETRON_INCLUDE_DIR}
    )

include_directories(${CMAKE_SOURCE_DIR}
                    ${ISMRMRD_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_PATH}
                    ${NUMPY_INCLUDE_DIRS}
                    )

set( core_header_files gadgetron_siemens_toolbox_core_export.h )
set( core_src_files )

set( application_header_files
    application/core_moco_averaging.h
    application/core_phase_sensitive_recon.h
    application/core_surface_coil_correction.h
    )

set( application_src_files 
    application/core_moco_averaging.cpp 
    application/core_phase_sensitive_recon.cpp 
    application/core_surface_coil_correction.cpp 
    )

set( schema_files )

source_group(schema FILES ${schema_files})

set(source_files ${core_header_files}
        ${core_src_files}
        ${application_header_files}
        ${application_src_files} )

source_group(application FILES  ${application_header_files} ${application_src_files})

add_library(gadgetron_siemens_toolbox_core SHARED ${source_files})

set_target_properties(gadgetron_siemens_toolbox_core PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_siemens_toolbox_core
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

target_include_directories(gadgetron_siemens_toolbox_core PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/application>
        )
install(TARGETS gadgetron_siemens_toolbox_core DESTINATION lib)
