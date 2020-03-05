
if (WIN32)
    add_definitions(-D__BUILD_SIEMENS_GADGETRON_SASHA_HC__)
endif ()

include_directories(
    ${ACE_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/toolboxes
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/application
    ${CMAKE_SOURCE_DIR}/toolboxes/cmr
    ${CMAKE_SOURCE_DIR}/toolboxes/cmr/application
    ${CMAKE_SOURCE_DIR}/gadgets/core
    ${CMAKE_SOURCE_DIR}/gadgets/sasha_hc
  )

INCLUDE_DIRECTORIES(${GADGETRON_INCLUDE_DIR})

IF (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
    ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

set( config_files
    config/SASHA-HC.xml 
    config/SASHA-HC_grappa.xml 
    config/SASHA-HC_grappa_moco.xml 
    config/SASHA-HC_grappa_moco_debug.xml 
    config/SASHA-HC_grappa_moco_istore.xml
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
     BucketToBufferHCGadget.h
     ReconGrappaSashaHCGadget.h
     SashaHCMultiRepGadget.h
     MoCoSashaHCGadget.h
     CmrParametricSashaT1T2MappingGadget.h )

set( sasha_hc_src_files 
     FlagSashaHCGadget.cpp
     CalcMolliTiGadget.cpp
     BucketToBufferHCGadget.cpp 
     ReconGrappaSashaHCGadget.cpp
     SashaHCMultiRepGadget.cpp
     MoCoSashaHCGadget.cpp
     CmrParametricSashaT1T2MappingGadget.cpp )

source_group(config     FILES      ${config_files})
source_group(schema     FILES      ${config_schema_files})
source_group(sasha_hc   FILES      ${sasha_hc_header_files} ${sasha_hc_src_files})

add_library(gadgetron_siemens_sasha_hc SHARED
            ${sasha_hc_header_files}
            ${sasha_hc_src_files}
            ${config_files}
            ${config_schema_files}
        )

set_target_properties(gadgetron_siemens_sasha_hc PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_siemens_sasha_hc
    gadgetron_siemens_toolbox_cmr 
    optimized gadgetron_cmr debug gadgetron_cmr${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpucore debug gadgetron_toolbox_cpucore${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpucore_math debug gadgetron_toolbox_cpucore_math${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpufft debug gadgetron_toolbox_cpufft${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpuklt debug gadgetron_toolbox_cpuklt${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_gadgettools debug gadgetron_toolbox_gadgettools${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_gadgetbase debug gadgetron_gadgetbase${CMAKE_DEBUG_SUFFIX}
    #optimized gadgetron_toolbox_cloudbus debug gadgetron_toolbox_cloudbus${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_log debug gadgetron_toolbox_log${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_mri_core debug gadgetron_toolbox_mri_core${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpudwt debug gadgetron_toolbox_cpudwt${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpuoperator debug gadgetron_toolbox_cpuoperator${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_mricore debug gadgetron_mricore${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpu_image debug gadgetron_toolbox_cpu_image${CMAKE_DEBUG_SUFFIX}
    ${Boost_LIBRARIES}
    ${ISMRMRD_LIBRARIES}
    ${FFTW3_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
    )

install (FILES  ${config_files}
                DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})

install (FILES  ${config_schema_files}
                DESTINATION ${GADGETRON_INSTALL_SCHEMA_PATH})

install(TARGETS gadgetron_siemens_sasha_hc DESTINATION lib)