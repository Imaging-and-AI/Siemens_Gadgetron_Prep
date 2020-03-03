#include "BucketToBufferHCGadget.h"

using BufferKey =  Gadgetron::BucketToBufferGadget::BufferKey;

namespace std {
    template<>
    struct less<BufferKey>{
        bool operator()(const BufferKey& idx1, const BufferKey& idx2) const {
            return std::tie(idx1.average,idx1.slice,idx1.contrast,idx1.phase,idx1.repetition,idx1.set,idx1.segment) <
                std::tie(idx2.average,idx2.slice,idx2.contrast,idx2.phase,idx2.repetition,idx2.set,idx2.segment);
        }
    };

    template<> struct equal_to<BufferKey>{
        bool operator()(const BufferKey& idx1, const BufferKey& idx2) const {
            return idx1.average == idx2.average
                   && idx1.slice == idx2.slice && idx1.contrast == idx2.contrast && idx1.phase == idx2.phase
                   && idx1.repetition == idx2.repetition && idx1.set == idx2.set && idx1.segment == idx2.segment;
        }
    };
}

namespace Gadgetron {

    namespace {

        IsmrmrdReconBit& getRBit(std::map<BufferKey, IsmrmrdReconData>& recon_data_buffers,
            const BufferKey& key, uint16_t espace) {

            // Look up the DataBuffered entry corresponding to this encoding space
            // create if needed and set the fields of view and matrix size
            if (recon_data_buffers[key].rbit_.size() < (espace + 1)) {
                recon_data_buffers[key].rbit_.resize(espace + 1);
            }

            return recon_data_buffers[key].rbit_[espace];
        }

    }

    BucketToBufferHCGadget::BucketToBufferHCGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : BucketToBufferGadget(context, props) {}

    // BucketToBufferHCGadget::BucketToBufferHCGadget()
    // {
    // }

    BucketToBufferHCGadget::~BucketToBufferHCGadget()
    {
        //The buckets array should be empty but just in case, let's make sure all the stuff is released.
    }

    void BucketToBufferHCGadget::process(Core::InputChannel<AcquisitionBucket>& input, Core::OutputChannel& out) {

        for (auto acq_bucket : input) {
            std::map<BufferKey, IsmrmrdReconData> recon_data_buffers;
            GDEBUG_STREAM("BUCKET_SIZE " << acq_bucket.data_.size() << " ESPACE " << acq_bucket.refstats_.size());
            // Iterate over the reference data of the bucket
            for (auto& acq : acq_bucket.ref_) {
                // Get a reference to the header for this acquisition

                const auto& acqhdr    = std::get<ISMRMRD::AcquisitionHeader>(acq);
                auto key              = getKey(acqhdr.idx);
                uint16_t espace       = acqhdr.encoding_space_ref;
                IsmrmrdReconBit& rbit = getRBit(recon_data_buffers, key, espace);
                if (!rbit.ref_) {
                    rbit.ref_ = makeDataBuffer(acqhdr, header.encoding[espace], acq_bucket.refstats_[espace], true);
                    rbit.ref_->sampling_ = createSamplingDescription(
                        header.encoding[espace], acq_bucket.refstats_[espace], acqhdr, true);
                }

                add_acquisition(*rbit.ref_, acq, header.encoding[espace], acq_bucket.refstats_[espace], true);

                // Stuff the data, header and trajectory into this data buffer
            }

            // Kelvin: Create a copy of reference data for second (high-contrast) encoding space
            for (auto& acq : acq_bucket.ref_) {
                // Get a reference to the header for this acquisition

                // NB: Not const anymore...
                auto& acqhdr    = std::get<ISMRMRD::AcquisitionHeader>(acq);

                acqhdr.encoding_space_ref = 1;

                auto key              = getKey(acqhdr.idx);
                uint16_t espace       = acqhdr.encoding_space_ref;
                IsmrmrdReconBit& rbit = getRBit(recon_data_buffers, key, espace);
                if (!rbit.ref_) {
                    rbit.ref_ = makeDataBuffer(acqhdr, header.encoding[espace], acq_bucket.refstats_[espace], true);
                    rbit.ref_->sampling_ = createSamplingDescription(
                        header.encoding[espace], acq_bucket.refstats_[espace], acqhdr, true);
                }

                add_acquisition(*rbit.ref_, acq, header.encoding[espace], acq_bucket.refstats_[espace], true);

                // Stuff the data, header and trajectory into this data buffer
            }

            // Kelvin: Find the limits of the high-contrast lines
            int16_t iHCLineMin = (int16_t)32767;
            int16_t iHCLineMax = 0;

            // Iterate over the data of the bucket
            for (auto& acq : acq_bucket.data_) {
                // Get a reference to the header for this acquisition

                const auto& acqhdr    = std::get<ISMRMRD::AcquisitionHeader>(acq);
                auto key              = getKey(acqhdr.idx);
                uint16_t espace       = acqhdr.encoding_space_ref;
                IsmrmrdReconBit& rbit = getRBit(recon_data_buffers, key, espace);
                if (rbit.data_.data_.empty()) {
                    rbit.data_ = makeDataBuffer(acqhdr, header.encoding[espace], acq_bucket.datastats_[espace], false);
                    rbit.data_.sampling_ = createSamplingDescription(
                        header.encoding[espace], acq_bucket.datastats_[espace], acqhdr, false);
                }

                add_acquisition(rbit.data_, acq, header.encoding[espace], acq_bucket.datastats_[espace], false);

                // Stuff the data, header and trajectory into this data buffer

                // 2016-10-03 Kelvin: Update limits for HC images
                if (acqhdr.idx.contrast)
                {
                    if ((int16_t)acqhdr.idx.kspace_encode_step_1 < iHCLineMin)
                    {
                        iHCLineMin = (int16_t)acqhdr.idx.kspace_encode_step_1;
                        //				GDEBUG_STREAM("Changed iHCLineMin: set: " << (int16_t)acqhdr.idx.set << " lin: " << (int16_t)acqhdr.idx.kspace_encode_step_1 << std::endl);
                    }

                    if ((int16_t)acqhdr.idx.kspace_encode_step_1 > iHCLineMax)
                    {
                        iHCLineMax = (int16_t)acqhdr.idx.kspace_encode_step_1;
                    }
                }
            }
            // Kelvin: Find the limits of the high-contrast lines
            GDEBUG_STREAM( "iHCLineMin: " << iHCLineMin << " iHCLineMax: " << iHCLineMax );

            // Kelvin: TODO: Should re-stuff the shared lines, but it was commented out before... WHY???

            // Send all the ReconData messages
            GDEBUG("End of bucket reached, sending out %d ReconData buffers\n", recon_data_buffers.size());

            for (auto& recon_data_buffer : recon_data_buffers) {
                if (acq_bucket.waveform_.empty())
                    out.push(recon_data_buffer.second);
                else
                    out.push(recon_data_buffer.second, acq_bucket.waveform_);
            }
        }
    }

    // GADGET_FACTORY_DECLARE( BucketToBufferHCGadget )
    GADGETRON_GADGET_EXPORT( BucketToBufferHCGadget )

}
