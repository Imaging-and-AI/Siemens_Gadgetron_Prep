#ifndef GT_MOCO_AVE_EXPORT_H_
#define GT_MOCO_AVE_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_SIEMENS_GADGETRON_MOCO_AVE__)
    #define EXPORTGTGADGET __declspec(dllexport)
#else
    #define EXPORTGTGADGET __declspec(dllimport)
#endif
#else
    #define EXPORTGTGADGET
#endif

#endif /* GT_MOCO_AVE_EXPORT_H_ */