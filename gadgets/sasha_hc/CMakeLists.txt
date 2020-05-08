
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
     FlagSashaHCGadget.h
     CalcMolliTiGadget.h
        ReconGrappaSashaHCGadget.h
     SashaHCMultiRepGadget.h
     MoCoSashaHCGadget.h
     CmrParametricSashaT1T2MappingGadget.h )

set( sasha_hc_src_files 
     FlagSashaHCGadget.cpp
     CalcMolliTiGadget.cpp
        ReconGrappaSashaHCGadget.cpp
     SashaHCMultiRepGadget.cpp
     MoCoSashaHCGadget.cpp
     CmrParametricSashaT1T2MappingGadget.cpp DuplicateReferenceDataGadget.cpp DuplicateReferenceDataGadget.h)

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
    Gadgetron::gadgetron_toolbox_cmr
    Gadgetron::gadgetron_cmr
    Gadgetron::gadgetron_toolbox_cpucore
    Gadgetron::gadgetron_toolbox_cpucore_math
    Gadgetron::gadgetron_toolbox_cpufft
    Gadgetron::gadgetron_toolbox_cpuklt
    Gadgetron::gadgetron_toolbox_log
    Gadgetron::gadgetron_toolbox_mri_core
    Gadgetron::gadgetron_toolbox_cpudwt
    Gadgetron::gadgetron_toolbox_cpuoperator
    Gadgetron::gadgetron_toolbox_cpu_image

    )

install (FILES  ${config_files}
                DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})

install (FILES  ${config_schema_files}
                DESTINATION ${GADGETRON_INSTALL_SCHEMA_PATH})

install(FILES ${python_files}
                DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH})

install(TARGETS gadgetron_siemens_sasha_hc DESTINATION lib)