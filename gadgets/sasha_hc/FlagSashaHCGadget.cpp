//FlagSashaHCGadget.cpp

#include "FlagSashaHCGadget.h"
#include <iomanip>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

using namespace Gadgetron;

FlagSashaHCGadget::FlagSashaHCGadget()
{
    new_sat_time_stamp_ = false;
    sat_time_stamp_     = 0;
    sat_time_stamp_set_ = -1;
}

int FlagSashaHCGadget::process_config( ACE_Message_Block* mb )
{
    GDEBUG_STREAM("FlagSashaHCGadget::process_config(...) starts ... ");
//    GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

    return GADGET_OK;
}

int FlagSashaHCGadget::process(
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

    // 2016-09-28 Kelvin: Mark high-contrast images with the contrasts dimension
    static int16_t iLastLine = -1;
    static int16_t iLastSet  = -1;

    // Reset iLastLine logic when a new image starts
    // Line reset to 0 is for ACS lines... FIXME!
    if (((int16_t)acqhdr.idx.set != iLastSet) || ((int16_t)acqhdr.idx.kspace_encode_step_1 == 0))
    {
        iLastSet = (int16_t)acqhdr.idx.set;
        iLastLine = -1;
    }

    // Assume that lines are always acquired in ascending order unless they are HC
    if ((int16_t)acqhdr.idx.kspace_encode_step_1 < iLastLine)
    {
        acqhdr.idx.contrast = 1;

        // 2017-02-15 Kelvin: Change the encoding space for the HC images
        acqhdr.encoding_space_ref = 1;
    }
    else if ( !(acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) || acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING)) ) {
        iLastLine = (int16_t)acqhdr.idx.kspace_encode_step_1;
    }

    if (acqhdr.idx.contrast)
    {
        GDEBUG_STREAM( "  set: "         <<                   (int16_t)acqhdr.idx.set
                    << " con: "         <<                   (int16_t)acqhdr.idx.contrast
                    << " e1: "          << std::setw( 3 ) << (int16_t)acqhdr.idx.kspace_encode_step_1
                    << " iLastLine: "   << std::setw( 3 ) << iLastLine
                    << " ts: "          << std::setw( 5 ) << (int16_t)acqhdr.user_int[4]
                    << " rfTime: "      << std::setw( 8 ) << (int32_t)acqhdr.user_int[7]
                    << " acquisition_time_stamp: "  <<       (uint32_t)acqhdr.acquisition_time_stamp
                    << " ref: "         <<                          ( acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) || acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING) )
                    << " *"); 
    }
    else {
        GDEBUG_STREAM( " set: "         <<                   (int16_t)acqhdr.idx.set
                    << " con: "         <<                   (int16_t)acqhdr.idx.contrast
                    << " e1: "          << std::setw( 3 ) << (int16_t)acqhdr.idx.kspace_encode_step_1
                    << " iLastLine: "   << std::setw( 3 ) << iLastLine
                    << " ts: "          << std::setw( 5 ) << (int16_t)acqhdr.user_int[4]
                    << " rfTime: "      << std::setw( 8 ) << (int32_t)acqhdr.user_int[7]
                    << " acquisition_time_stamp: "  <<       (uint32_t)acqhdr.acquisition_time_stamp
                    << " ref: "         <<                          ( acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) || acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING) ) ); 
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

GADGET_FACTORY_DECLARE( FlagSashaHCGadget )