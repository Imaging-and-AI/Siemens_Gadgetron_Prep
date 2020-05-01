//
// Created by david on 29/04/2020.
//

#include "DuplicateReferenceDataGadget.h"
Gadgetron::AcquisitionBucket Gadgetron::DuplicateReferenceDataGadget::process_function(
    Gadgetron::AcquisitionBucket bucket) const {

    //Copy the reference data
    auto new_reference_data = std::vector<Core::Acquisition>{};
    new_reference_data.reserve(bucket.ref_.size());
    for (auto& [header,data,traj] : bucket.ref_){
        if (header.encoding_space_ref == encoding_space_from) new_reference_data.emplace_back(header,data,traj);
    }

    GDEBUG("Copying %d reference lines. Total is %d\n",new_reference_data.size(),bucket.ref_.size());

    //Adjust the encoding_space_ref
    for (auto& [header,data,traj] : new_reference_data){
         header.encoding_space_ref = encoding_space_to;
    }
    //Append to the original reference data
    bucket.ref_.insert(bucket.ref_.end(),new_reference_data.begin(),new_reference_data.end());

    bucket.refstats_.resize(std::max<size_t>(bucket.refstats_.size(),encoding_space_to+1));
    bucket.refstats_[encoding_space_to] = bucket.refstats_[0];

    return std::move(bucket);

}
Gadgetron::DuplicateReferenceDataGadget::DuplicateReferenceDataGadget(
    const Core::Context& context, const Core::GadgetProperties& props)
    : Core::PureGadget<AcquisitionBucket,AcquisitionBucket>(props) { }

namespace Gadgetron {
    GADGETRON_GADGET_EXPORT(DuplicateReferenceDataGadget)
}
