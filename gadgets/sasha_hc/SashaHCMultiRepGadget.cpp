
#include "SashaHCMultiRepGadget.h"
#include <iomanip>
#include <sstream>

namespace Gadgetron { 

    SashaHCMultiRepGadget::SashaHCMultiRepGadget() : BaseClass()
    {
    }

    SashaHCMultiRepGadget::~SashaHCMultiRepGadget()
    {
    }

    int SashaHCMultiRepGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
        max_rep_ = e_limits.repetition ? e_limits.repetition->maximum : 0;

        GDEBUG_STREAM("Maximal repetition is " << max_rep_);

        buf_lc_.resize(max_rep_ + 1);
        buf_hc_.resize(max_rep_ + 1);

        return GADGET_OK;
    }

    int SashaHCMultiRepGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("SashaHCMultiRepGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "SashaHCMultiRepGadget::process(...) starts ... ");

        // -------------------------------------------------------------

        process_called_times_++;

        // -------------------------------------------------------------

        IsmrmrdImageArray* data = m1->getObjectPtr();

        // print out data info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> SashaHCMultiRepGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            data->data_.print(os);
            GDEBUG_STREAM(os.str());
        }

        // -------------------------------------------------------------

        // some images do not need mapping
        if (data->meta_[0].length(skip_processing_meta_field.value().c_str())>0)
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("SashaHCMultiRepGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        // -------------------------------------------------------------

        size_t encoding = (size_t)data->meta_[0].as_long("encoding", 0);

        std::string dataRole = std::string(data->meta_[0].as_str(GADGETRON_DATA_ROLE));

        size_t RO = data->data_.get_size(0);
        size_t E1 = data->data_.get_size(1);
        size_t E2 = data->data_.get_size(2);
        size_t CHA = data->data_.get_size(3);
        size_t N = data->data_.get_size(4);
        size_t S = data->data_.get_size(5);
        size_t SLC = data->data_.get_size(6);

        std::stringstream os;
        os << "dataRole_" << dataRole << "_encoding_" << encoding << "_processing_" << process_called_times_;
        std::string str = os.str();

        // -------------------------------------------------------------

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(data->data_, debug_folder_full_path_ + "data" + str);
        }

        // check incoming data rep
        int32_t con = m1->getObjectPtr()->headers_(0, 0, 0).contrast;
        int32_t rep = m1->getObjectPtr()->headers_(0, 0, 0).repetition;

        //if(max_rep_==0)
        //{
        //    if (this->next()->putq(m1) == -1)
        //    {
        //        GERROR("SashaHCMultiRepGadget::process, passing data on to next gadget");
        //        return GADGET_OK;
        //    }
        //}

        if (con == 0)
            buf_lc_[rep] = *m1->getObjectPtr();
        else
            buf_hc_[rep] = *m1->getObjectPtr();

        if (rep == max_rep_ && con==1)
        {
            size_t REP = max_rep_ + 1;

            // combine and send out data
            Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm1 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();

            IsmrmrdImageArray* pRes = cm1->getObjectPtr();

            pRes->data_.create(RO, E1, E2, CHA, N*REP, 2*S, SLC);
            pRes->headers_.create(N*REP, 2*S, SLC);
            pRes->meta_.resize(N*REP* 2 * S*SLC);

            size_t slc, n, rep;
            for (slc = 0; slc<SLC; slc++)
            {
                for (rep = 0; rep < REP; rep++)
                {
                    memcpy(&pRes->data_(0, 0, 0, 0, rep*N, 0, slc), &buf_lc_[rep].data_(0, 0, 0, 0, 0, 0, slc), sizeof(std::complex<float>)*RO*E1*E2*CHA*N);
                    memcpy(&pRes->data_(0, 0, 0, 0, rep*N, 1, slc), &buf_hc_[rep].data_(0, 0, 0, 0, 0, 0, slc), sizeof(std::complex<float>)*RO*E1*E2*CHA*N);

                    for (n = 0; n < N; n++)
                    {
                        pRes->headers_(n + rep*N, 0, slc) = buf_lc_[rep].headers_(n, 0, slc);
                        pRes->meta_[n + rep*N + 0*N*REP + slc*N*S*REP] = buf_lc_[rep].meta_[n+0*N+slc*N*S];

                        pRes->headers_(n + rep*N, 1, slc) = buf_hc_[rep].headers_(n, 0, slc);
                        pRes->meta_[n + rep*N + 1 * N*REP + slc*N*S*REP] = buf_hc_[rep].meta_[n + 0 * N + slc*N*S];
                    }
                }
            }

            if (this->next()->putq(cm1) == -1)
            {
                GERROR("SashaHCMultiRepGadget::process, passing data on to next gadget");
                return GADGET_FAIL;
            }
        }

        m1->release();

        // -------------------------------------------------------------

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    int SashaHCMultiRepGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "SashaHCMultiRepGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }


    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(SashaHCMultiRepGadget)

}
