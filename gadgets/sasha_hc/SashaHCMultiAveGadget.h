/** \file   SashaHCMultiAveGadget.h
    \brief  This is the class gadget to handle multi-rep sasha for high contrast use case.
    \author Hui Xue
*/

#pragma once

#include "sashahclib_export.h"
#include "GenericReconBase.h"

namespace Gadgetron {

    class EXPORTSASHAHC SashaHCMultiAveGadget : public Gadgetron::GenericReconImageBase
    {
    public:
        GADGET_DECLARE(SashaHCMultiAveGadget);

        typedef Gadgetron::GenericReconImageBase BaseClass;

        SashaHCMultiAveGadget();
        ~SashaHCMultiAveGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the mapping
        /// ------------------------------------------------------------------------------------

        GADGET_PROPERTY(skip_processing_meta_field, std::string, "If this meta field exists, pass the incoming image array to next gadget without processing", GADGETRON_SKIP_PROCESSING_AFTER_RECON);

        // ------------------------------------------------------------------------------------

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        int max_ave_;

        std::vector<IsmrmrdImageArray> buf_lc_;
        std::vector<IsmrmrdImageArray> buf_hc_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

        // close call
        int close(unsigned long flags);
    };
}
