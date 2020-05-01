#ifndef RECONGRAPPASASHAHCGRAPPA_H
#define RECONGRAPPASASHAHCGRAPPA_H

#include "sashahclib_export.h"
#include <gadgetron/GenericReconCartesianGrappaGadget.h>

namespace Gadgetron {
    class EXPORTSASHAHC ReconGrappaSashaHCGadget :
        public GenericReconCartesianGrappaGadget
    {
    public:
        GADGET_DECLARE( ReconGrappaSashaHCGadget );

        ReconGrappaSashaHCGadget();
        virtual ~ReconGrappaSashaHCGadget();

    protected:
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);
        virtual void perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e);

        // Storage of GRAPPA-ed k-space data from previous encoding space
        hoNDArray< std::complex<float> > unaliased_uncombined_enc0_;
    };
}
#endif //RECONGRAPPASASHAHCGRAPPA_H
