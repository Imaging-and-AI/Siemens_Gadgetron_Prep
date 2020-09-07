//ThresholdGadget.h

#ifndef FLAGSASHAHCGADGET_H
#define FLAGSASHAHCGADGET_H

#include "sashahclib_export.h"
#include <gadgetron/Gadget.h>
#include <gadgetron/GadgetMRIHeaders.h>
#include <gadgetron/hoNDArray.h>
#include <complex>
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{

    class EXPORTSASHAHC FlagSashaHCGadget :
        public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE( FlagSashaHCGadget )

        FlagSashaHCGadget();

        GADGET_PROPERTY(verbose, bool, "verbose mode", false);

    protected:
        virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 );

        virtual int process_config( ACE_Message_Block* mb );

        //  float threshold_level_;
        bool     new_sat_time_stamp_;
        uint32_t sat_time_stamp_;
        int16_t  sat_time_stamp_set_;

    };

}
#endif //FLAGSASHAHCGADGET_H