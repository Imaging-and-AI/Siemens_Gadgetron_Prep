//
#pragma once
#include <gadgetron/PureGadget.h>
#include <gadgetron/mri_core_acquisition_bucket.h>
#include <gadgetron/mri_core_data.h>
namespace Gadgetron {
class DuplicateReferenceDataGadget
    : public Core::PureGadget<AcquisitionBucket, AcquisitionBucket> {
public:
    DuplicateReferenceDataGadget(const Core::Context& context, const Core::GadgetProperties& props);
  AcquisitionBucket process_function(AcquisitionBucket args) const override;

  NODE_PROPERTY(encoding_space_to,uint16_t,"Encoding space to copy data to",1);
  NODE_PROPERTY(encoding_space_from,uint16_t,"Encoding space to copy data to",0);
};
} // namespace Gadgetron
