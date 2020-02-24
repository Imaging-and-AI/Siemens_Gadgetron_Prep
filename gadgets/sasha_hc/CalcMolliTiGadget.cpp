//CalcMolliTiGadget.cpp

#include "CalcMolliTiGadget.h"
#include <iomanip>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

using namespace Gadgetron;

CalcMolliTiGadget::CalcMolliTiGadget()
{
    ref_last_ksp_timestamp_  = 0;
    curr_last_ksp_timestamp_ = 0;
    curr_look_locker_set_     = -1;
}

int CalcMolliTiGadget::process_config( ACE_Message_Block* mb )
{
    return GADGET_OK;
}

int CalcMolliTiGadget::process(
    GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 )
{
    ISMRMRD::AcquisitionHeader & acqhdr = *m1->getObjectPtr();

    // Ignore noise and ref lines
    if ( !acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT) && 
        !(acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) || acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING)) )
    {
        if (acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
        {
            curr_last_ksp_timestamp_ = acqhdr.acquisition_time_stamp;

            if (acqhdr.user_int[5] != curr_look_locker_set_)
            {
                // First line of first image in Look-Locker set
                ref_last_ksp_timestamp_  = acqhdr.acquisition_time_stamp;
                curr_look_locker_set_ = acqhdr.user_int[5];
            }

            // Store the TI offset from the first image in LL set
            acqhdr.user_int[7] = (curr_last_ksp_timestamp_ - ref_last_ksp_timestamp_)*2500;

            // Add the min TI of the current Look-Locker set to get actual TI
            acqhdr.user_int[7] += acqhdr.user_int[4]*1000;

            GDEBUG_STREAM( " set: "            <<                    (int16_t)acqhdr.idx.set
                        << " e1: "             << std::setw( 3 ) <<  (int16_t)acqhdr.idx.kspace_encode_step_1
                        << " ti: "             << std::setw( 3 ) <<  (int16_t)acqhdr.user_int[4]
                        << " ll: "             << std::setw( 1 ) <<  (int16_t)acqhdr.user_int[5]
                        << " tieff: "          << std::setw( 8 ) << (uint32_t)acqhdr.user_int[7]
                        << " acq_time_stamp: " << std::setw( 8 ) << (uint32_t)acqhdr.acquisition_time_stamp
                        << " ref: "            << std::setw( 8 ) << (uint32_t)ref_last_ksp_timestamp_
                        );
        }
    }

    //Now pass on image
    if (this->next()->putq( m1 ) < 0) {
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE( CalcMolliTiGadget )