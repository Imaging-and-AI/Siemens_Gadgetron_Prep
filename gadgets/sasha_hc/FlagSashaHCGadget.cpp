//FlagSashaHCGadget.cpp

#include "FlagSashaHCGadget.h"
#include <iomanip>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

using namespace Gadgetron;

int FlagSashaHCGadget::process_config( ACE_Message_Block* mb )
{
    return GADGET_OK;
}

int FlagSashaHCGadget::process(
    GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 )
{

    bool is_noise = ISMRMRD::FlagBit( ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT ).isSet( m1->getObjectPtr()->flags );

    if (!is_noise)
    {

        // 2016-09-28 Kelvin: Mark high-contrast images with the contrasts dimension
        ISMRMRD::AcquisitionHeader & acqhdr = *m1->getObjectPtr();
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
            GDEBUG_STREAM( " set: "         <<                   (int16_t)acqhdr.idx.set
                        << " con: "         <<                   (int16_t)acqhdr.idx.contrast
                        << " e1: "          << std::setw( 3 ) << (int16_t)acqhdr.idx.kspace_encode_step_1
                        << " iLastLine: "   << std::setw( 3 ) << iLastLine
                        << " ts: "          <<                   (int16_t)acqhdr.user_int[4]
                        << " ref: "         <<                          ( acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) || acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING) )
                        << " *"); 
        }
        else {
            GDEBUG_STREAM( " set: "         <<                   (int16_t)acqhdr.idx.set
                        << " con: "         <<                   (int16_t)acqhdr.idx.contrast
                        << " e1: "          << std::setw( 3 ) << (int16_t)acqhdr.idx.kspace_encode_step_1
                        << " iLastLine: "   << std::setw( 3 ) << iLastLine
                        << " ts: "          <<                   (int16_t)acqhdr.user_int[4]
                        << " ref: "         <<                          ( acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) || acqhdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING) ) ); 
        }

        // 2016-10-04 Kelvin: Discard navigator lines
        if (acqhdr.isFlagSet( ISMRMRD::ISMRMRD_ACQ_IS_RTFEEDBACK_DATA ))
        {
            return GADGET_OK;
        }
    }

    //Now pass on image
    if (this->next()->putq( m1 ) < 0) {
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE( FlagSashaHCGadget )