/** \file       gadgetron_siemens_toolbox_core_export.h
    \brief      Implement windows export/import for core toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_SIEMENS_TOOLBOX_CORE__) || defined (gadgetron_siemens_toolbox_core_EXPORTS)
        #define EXPORTGTTOOLBOXCORE __declspec(dllexport)
    #else
        #define EXPORTGTTOOLBOXCORE __declspec(dllimport)
    #endif
#else
    #define EXPORTGTTOOLBOXCORE
#endif
