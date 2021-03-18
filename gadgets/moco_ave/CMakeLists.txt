
if (WIN32)
    add_definitions(-D__BUILD_SIEMENS_GADGETRON_MOCO_AVE__)
endif ()

IF (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
    ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

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
                    ${ARMADILLO_LIBRARIES} 
                    )

set( config_files
    config/Siemens_Gadgetron_LGE_MOCO_AVE.xml 
    config/Siemens_Gadgetron_LGE_MOCO_AVE_dstore.xml 
    config/Siemens_Gadgetron_LGE_MOCO_AVE_istore.xml 
    config/Siemens_Gadgetron_WB_LGE_MOCO_AVE.xml 
    config/Siemens_Gadgetron_DB_LGE_MOCO_AVE_istore.xml
    )

set( config_schema_files )

set( gt_header_files 
     gadgetron_siemens_moco_ave_export.h 
     GtImageMoCoAveGadget.h 
     GtImageMoCoAvePSIRGadget.h )

set( gt_src_files 
     GtImageMoCoAveGadget.cpp 
     GtImageMoCoAvePSIRGadget.cpp )

source_group(config     FILES      ${config_files})
source_group(schema     FILES      ${config_schema_files})
source_group(gt         FILES      ${gt_header_files} ${gt_src_files})

add_library(gadgetron_siemens_moco_ave SHARED
            ${gt_header_files}
            ${gt_src_files}
            ${config_files}
        )


target_link_libraries(gadgetron_siemens_moco_ave
    gadgetron_siemens_toolbox_core 
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpucore_math
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpuklt
    gadgetron_toolbox_log
    gadgetron_toolbox_mri_core
    gadgetron_toolbox_cpudwt
    gadgetron_toolbox_cpuoperator
    gadgetron_mricore
    gadgetron_toolbox_cpu_image
    gadgetron_toolbox_fatwater
    gadgetron_toolbox_denoise
    gadgetron_toolbox_image_analyze_io

    gadgetron_core
    gadgetron_core_writers
    gadgetron_core_readers
    gadgetron_core_parallel
    gadgetron_core_distributed

    gadgetron_toolbox_image_analyze_io
    gadgetron_toolbox_hostutils

    ${ISMRMRD_LIBRARIES}
    Boost::system
    Boost::filesystem
    Boost::program_options
    )

install (FILES  ${config_files}
                DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})

install(TARGETS gadgetron_siemens_moco_ave DESTINATION lib)