//ThresholdGadget.h

#ifndef FLAGSASHAGADGET_H
#define FLAGSASHAGADGET_H

#include "sashahclib_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include <complex>
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{

    class EXPORTSASHAHC FlagSashaGadget :
        public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE( FlagSashaGadget )

        FlagSashaGadget();

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
#endif //FLAGSASHAGADGET_H