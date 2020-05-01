//ThresholdGadget.h

#ifndef CALCMOLLITIGADGET_H
#define CALCMOLLITIGADGET_H

#include "sashahclib_export.h"
#include <gadgetron/Gadget.h>
#include <gadgetron/GadgetMRIHeaders.h>
#include <gadgetron/hoNDArray.h>
#include <complex>
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{

    class EXPORTSASHAHC CalcMolliTiGadget :
        public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE( CalcMolliTiGadget )

        CalcMolliTiGadget();

    protected:
        virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 );

        virtual int process_config( ACE_Message_Block* mb );

        uint32_t ref_last_ksp_timestamp_;
        uint32_t curr_last_ksp_timestamp_;
        int16_t  curr_look_locker_set_;
    };

}
#endif //CALCMOLLITIGADGET_H