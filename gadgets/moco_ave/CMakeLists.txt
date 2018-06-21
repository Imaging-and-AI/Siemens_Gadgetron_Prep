
if (WIN32)
    add_definitions(-D__BUILD_SIEMENS_GADGETRON_MOCO_AVE__)
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
    ${CMAKE_SOURCE_DIR}/gadgets
    ${CMAKE_SOURCE_DIR}/gadgets/moco_ave
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
    config/Siemens_Gadgetron_LGE_MOCO_AVE.xml 
    config/Siemens_Gadgetron_LGE_MOCO_AVE_dstore.xml 
    config/Siemens_Gadgetron_LGE_MOCO_AVE_istore.xml 
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

set_target_properties(gadgetron_siemens_moco_ave PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_siemens_moco_ave
    gadgetron_siemens_toolbox_core 
    optimized gadgetron_toolbox_cpucore debug gadgetron_toolbox_cpucore${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpucore_math debug gadgetron_toolbox_cpucore_math${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpufft debug gadgetron_toolbox_cpufft${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpuklt debug gadgetron_toolbox_cpuklt${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_gadgettools debug gadgetron_toolbox_gadgettools${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_gadgetbase debug gadgetron_gadgetbase${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cloudbus debug gadgetron_toolbox_cloudbus${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_log debug gadgetron_toolbox_log${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_mri_core debug gadgetron_toolbox_mri_core${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpudwt debug gadgetron_toolbox_cpudwt${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpuoperator debug gadgetron_toolbox_cpuoperator${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_mricore debug gadgetron_mricore${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_cpu_image debug gadgetron_toolbox_cpu_image${CMAKE_DEBUG_SUFFIX}
    optimized gadgetron_toolbox_image_analyze_io debug gadgetron_toolbox_image_analyze_io${CMAKE_DEBUG_SUFFIX}
    ${ISMRMRD_LIBRARIES}
    ${FFTW3_LIBRARIES} 
    ${Boost_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
    )

install (FILES  ${config_files}
                DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})

install(TARGETS gadgetron_siemens_moco_ave DESTINATION lib)