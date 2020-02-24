/** \file       gadgetron_siemens_toolbox_cmr_export.h
    \brief      Implement windows export/import for cmr toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_SIEMENS_TOOLBOX_CMR__) || defined (gadgetron_siemens_toolbox_cmr_EXPORTS)
        #define EXPORTGTTOOLBOXCMR __declspec(dllexport)
    #else
        #define EXPORTGTTOOLBOXCMR __declspec(dllimport)
    #endif
#else
    #define EXPORTGTTOOLBOXCMR
#endif
