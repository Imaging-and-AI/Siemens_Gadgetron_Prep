
if (WIN32)
    add_definitions(-D__BUILD_SIEMENS_GADGETRON_SASHA_HC__)
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
                    )

set(python_files
        MocoAntsHC.py
        )

set( config_files
    config/SASHA-HC.xml 
    config/SASHA-HC_grappa.xml 
    config/SASHA-HC_grappa_moco.xml
    config/SASHA-HC_grappa_nomoco.xml
    config/SASHA-HC_grappa_base.xml
    config/SASHA-HC_grappa_ants.xml
    config/SASHA-HC_ants.xml
    config/SASHA-HC_grappa_moco_debug.xml
    config/SASHA-HC_grappa_moco_istore.xml
    config/SASHA-HC_grappa_moco_dstore.xml
    )

set( config_schema_files
    ../../schema/IsmrmrdParameterMap_Siemens_SASHA_HC.xsl
    ../../schema/IsmrmrdParameterMap_Siemens_SashaT1T2.xsl
    ../../schema/IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xml
    ../../schema/IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xsl
    )

set( sasha_hc_header_files 
     sashahclib_export.h 
     CalcMolliTiGadget.h
     CmrParametricSashaT1T2MappingGadget.h
     DuplicateReferenceDataGadget.h
     FlagSashaHCGadget.h
     FlagSashaGadget.h
     MoCoSashaHCGadget.h
     MoCoSashaGadget.h
     ReconGrappaSashaHCGadget.h
     SashaHCMultiAveGadget.h
     SashaHCMultiRepGadget.h 
     )

set( sasha_hc_src_files 
    CalcMolliTiGadget.cpp
    CmrParametricSashaT1T2MappingGadget.cpp
    DuplicateReferenceDataGadget.cpp
    FlagSashaHCGadget.cpp
    FlagSashaGadget.cpp
    MoCoSashaHCGadget.cpp
    MoCoSashaGadget.cpp
    ReconGrappaSashaHCGadget.cpp
    SashaHCMultiAveGadget.cpp
    SashaHCMultiRepGadget.cpp 
    )

source_group(config     FILES      ${config_files})
source_group(schema     FILES      ${config_schema_files})
source_group(sasha_hc   FILES      ${sasha_hc_header_files} ${sasha_hc_src_files})

add_library(gadgetron_siemens_sasha_hc SHARED
            ${sasha_hc_header_files}
            ${sasha_hc_src_files}
            ${config_files}
            ${config_schema_files}
        )


target_link_libraries(gadgetron_siemens_sasha_hc
    gadgetron_siemens_toolbox_cmr 
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpucore_math
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpuklt
    gadgetron_toolbox_log
    gadgetron_toolbox_mri_core
    gadgetron_toolbox_cpudwt
    gadgetron_toolbox_cpuoperator
    gadgetron_toolbox_cmr
    gadgetron_mricore
    gadgetron_toolbox_cpu_image
    gadgetron_toolbox_fatwater
    gadgetron_toolbox_denoise
    gadgetron_toolbox_image_analyze_io
    gadgetron_cmr

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

install (FILES  ${config_schema_files}
                DESTINATION ${GADGETRON_INSTALL_SCHEMA_PATH})

install(FILES ${python_files}
                DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH})

install(TARGETS gadgetron_siemens_sasha_hc DESTINATION lib)