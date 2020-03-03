#ifndef BUCKETTOBUFFERHC_H
#define BUCKETTOBUFFERHC_H

#include "sashahclib_export.h"
#include "Gadget.h"
#include "BucketToBufferGadget.h"

namespace Gadgetron {

    class EXPORTSASHAHC BucketToBufferHCGadget :
        public BucketToBufferGadget
    {
    public:
        GADGET_DECLARE( BucketToBufferHCGadget );

        // BucketToBufferHCGadget();
        BucketToBufferHCGadget(const Core::Context& context, const Core::GadgetProperties& props);
        virtual ~BucketToBufferHCGadget();

    protected:
        void process(Core::InputChannel<AcquisitionBucket>& input, Core::OutputChannel& out) override;
    };
}
#endif //BUCKETTOBUFFERHC_H
