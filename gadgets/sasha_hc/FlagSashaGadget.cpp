//FlagSashaGadget.cpp

#include "FlagSashaGadget.h"
#include <iomanip>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

using namespace Gadgetron;

FlagSashaGadget::FlagSashaGadget()
{
    new_sat_time_stamp_ = false;
    sat_time_stamp_     = 0;
    sat_time_stamp_set_ = -1;
}

int FlagSashaGadget::process_config( ACE_Message_Block* mb )
{
    GDEBUG_STREAM("FlagSashaGadget::process_config(...) starts ... ");
//    GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

    return GADGET_OK;
}

int FlagSashaGadget::process(
    GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 )
{
    ISMRMRD::AcquisitionHeader & acqhdr = *m1->getObjectPtr();

    // Process dummy ADC for multi-HB sat recovery timing
    if ( acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT) && (acqhdr.user_int[6] != 0) )
    {
        new_sat_time_stamp_ = true;
        sat_time_stamp_     = acqhdr.acquisition_time_stamp;

        GDEBUG_STREAM("  Dummy ADC for multi-HB timing" << "                                   "
                   << "acquisition_time_stamp: "       <<  sat_time_stamp_   ); 

        // Do not pass this dummy ADC down the chain as it may interfere with other noise processing
        return GADGET_OK;
    }

    if (new_sat_time_stamp_)
    {
        // First line after the dummy ADC indicates the image set index that it applies to
        sat_time_stamp_set_ = (int16_t)acqhdr.idx.set;
        new_sat_time_stamp_ = false;
    }

    if ( (int16_t)acqhdr.idx.set == sat_time_stamp_set_ )
    {
        // Store the "timeSinceLastRF" in user_int[7], to be retrieved by CmrParametricSashaT1T2MappingGadget
        acqhdr.user_int[7] = (acqhdr.acquisition_time_stamp - sat_time_stamp_)*2500;
    }
    else
    {
        // Explicitly discard the stored value as it may be inaccurate
        acqhdr.user_int[7] = 0;
    }

    // 2016-10-04 Kelvin: Discard navigator lines
    if (acqhdr.isFlagSet( ISMRMRD::ISMRMRD_ACQ_IS_RTFEEDBACK_DATA ))
    {
        return GADGET_OK;
    }

    //Now pass on image
    if (this->next()->putq( m1 ) < 0) {
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE( FlagSashaGadget )