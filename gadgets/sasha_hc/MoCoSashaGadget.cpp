
#include "MoCoSashaGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "cmr_motion_correction.h"

namespace Gadgetron {

    MoCoSashaGadget::MoCoSashaGadget() : BaseClass()
    {
        this->do_moco_ = false;
    }

    MoCoSashaGadget::~MoCoSashaGadget()
    {
    }

    int MoCoSashaGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(),h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        num_encoding_spaces_ = 1;//h.encoding.size();
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << num_encoding_spaces_);

        iters_.resize(4);
        iters_[0] = 64;
        iters_[1] = 100;
        iters_[2] = 100;
        iters_[3] = 100;

        bool has_moco_in_proto = false;

        if (h.userParameters->userParameterLong.size() > 0)
        {
            std::vector<ISMRMRD::UserParameterLong>::const_iterator iter = h.userParameters->userParameterLong.begin();

            for (; iter != h.userParameters->userParameterLong.end(); iter++)
            {
                std::string usrParaName = iter->name;
                long        usrParaValue = iter->value;

                if (usrParaName == "MotionCorrection")
                {
                    this->do_moco_ = (usrParaValue>1);
                    has_moco_in_proto = true;
                    GDEBUG_STREAM("MoCoSashaGadget, found MotionCorrection in protocol : " << usrParaValue << " - do_moco is " << this->do_moco_);
                }
            }
        }

        if (!has_moco_in_proto)
        {
            this->do_moco_ = true;
            GDEBUG_STREAM("MoCoSashaGadget, cannot find MotionCorrection in protocol : " << this->do_moco_);
        }

        return GADGET_OK;
    }

    int MoCoSashaGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("MoCoSashaGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "MoCoSashaGadget::process(...) starts ... ");

        process_called_times_++;

        std::stringstream os;
        os << "called_" << process_called_times_;
        std::string str = os.str();

        IsmrmrdImageArray recon_res_ = *(m1->getObjectPtr());

        if (send_ori.value())
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("MoCoSashaGadget::process, passing m1 on to next gadget");
                return GADGET_FAIL;
            }
        }

        size_t encoding = (size_t)recon_res_.meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding < num_encoding_spaces_, GADGET_FAIL);

        std::string imageInfo;

        size_t RO = recon_res_.data_.get_size(0);
        size_t E1 = recon_res_.data_.get_size(1);
        size_t E2 = recon_res_.data_.get_size(2);
        size_t CHA = recon_res_.data_.get_size(3);
        size_t N = recon_res_.data_.get_size(4); // set
        size_t S = recon_res_.data_.get_size(5); // average
        size_t SLC = recon_res_.data_.get_size(6);

        if (!debug_folder_full_path_.empty())
        {
            hoNDArray<T> data = recon_res_.data_;
            data.squeeze();
            gt_exporter_.export_array_complex(data, debug_folder_full_path_ + "MoCoSasha_data_" + str);
        }

        size_t ro, e1, n, s, slc;

        bool warp_input = false;

        hoNDArray<T> moco;
        moco.create(RO, E1, E2, CHA, N, S, SLC);

        for (slc=0; slc<SLC; slc++)
        {
            std::stringstream os;
            os << "slice_" << slc;
            std::string str_slc = os.str();

            hoNDArray<real_value_type> mag, data_real, data_imag;
            mag.create(RO, E1, N*S);
            data_real.create(RO, E1, N*S);
            data_imag.create(RO, E1, N*S);

            for (s=0; s<S; s++)
            {
                for (n=0; n<N; n++)
                {
                    for (e1=0; e1<E1; e1++)
                    {
                        for (ro=0; ro<RO; ro++)
                        {
                            mag(ro, e1, n+s*N) = std::abs(recon_res_.data_(ro, e1, 0, 0, n, s, slc));
                            data_real(ro, e1, n+s*N) = recon_res_.data_(ro, e1, 0, 0, n, s, slc).real();
                            data_imag(ro, e1, n+s*N) = recon_res_.data_(ro, e1, 0, 0, n, s, slc).imag();
                        }
                    }
                }
            }

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(mag, debug_folder_full_path_ + "MoCoSasha_data_mag_" + str);
                gt_exporter_.export_array(data_real, debug_folder_full_path_ + "MoCoSasha_data_real_" + str);
                gt_exporter_.export_array(data_imag, debug_folder_full_path_ + "MoCoSasha_data_imag_" + str);
            }

            size_t key_frame(0);
            Gadgetron::find_key_frame_SSD_2DT(mag, key_frame);

            // perform moco (on diff images)
            std::string moco_str = "MoCoSashaGadget, perform moco for " + str_slc;
            if (perform_timing.value()) { gt_timer_local_.start(moco_str.c_str()); }
            Gadgetron::hoImageRegContainer2DRegistration<ImageType, ImageType, double> reg;

            Gadgetron::perform_moco_fixed_key_frame_2DT(mag, key_frame, (real_value_type)(10.0), iters_, bidirectional_moco.value(), warp_input, reg);

            if (perform_timing.value()) { gt_timer_local_.stop(); }

            // apply the deformation (on all images)
            moco_str = "MoCoSashaGadget, apply deformation field for " + str_slc;
            if (perform_timing.value()) { gt_timer_local_.start(moco_str.c_str()); }

            hoNDArray<double> dx, dy;
            reg.deformation_field_[0].to_NDArray(0, dx);
            reg.deformation_field_[1].to_NDArray(0, dy);

            hoNDArray<real_value_type> data_real_moco, data_imag_moco, data_mag_moco;
            Gadgetron::apply_deformation_field(data_real, dx, dy, data_real_moco);
            Gadgetron::apply_deformation_field(data_imag, dx, dy, data_imag_moco);
            Gadgetron::apply_deformation_field(mag, dx, dy, data_mag_moco);

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "called_" << process_called_times_ << "_slc" << slc;
                std::string str = os.str();

                gt_exporter_.export_array(data_real_moco, debug_folder_full_path_ + "MoCoSasha_data_real_moco_" + str);
                gt_exporter_.export_array(data_imag_moco, debug_folder_full_path_ + "MoCoSasha_data_imag_moco_" + str);
                gt_exporter_.export_array(data_mag_moco, debug_folder_full_path_ + "MoCoSasha_data_mag_moco_" + str);
            }

            hoNDArray<T> data_moco;
            data_moco.create(RO, E1, N*S);
            Gadgetron::real_imag_to_complex(data_real_moco, data_imag_moco, data_moco);

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "called_" << process_called_times_ << "_slc" << slc;
                std::string str = os.str();

                gt_exporter_.export_array_complex(data_moco, debug_folder_full_path_ + "MoCoSasha_data_moco_" + str);
            }

            for (s=0; s<S; s++)
            {
                for (n=0; n<N; n++)
                {
                    for (e1=0; e1<E1; e1++)
                    {
                        for (ro=0; ro<RO; ro++)
                        {
                            moco(ro, e1, 0, 0, n, s, slc) = data_moco(ro, e1, n+s*N);
                        }
                    }
                }
            }
        }

        if (!debug_folder_full_path_.empty())
        {
            std::stringstream os;
            os << "called_" << process_called_times_ << "_slc" << slc;
            std::string str = os.str();

            hoNDArray<T> data;
            data = moco; 
            data.squeeze();
            gt_exporter_.export_array_complex(data, debug_folder_full_path_ + "MoCoSasha_moco_res_" + str);
        }

        GDEBUG_CONDITION_STREAM(verbose.value(), "MoCoSashaGadget::process(...) ends ... ");

        // ----------------------------------------------------
        // format data to be passed on
        // ----------------------------------------------------
        Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_moco      = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();

        if (send_moco.value())
        {
            m1_sasha_moco->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, S, SLC);
            m1_sasha_moco->getObjectPtr()->headers_.create(N, S, SLC);
            m1_sasha_moco->getObjectPtr()->meta_.resize(N*S*SLC);
        }

        if(this->disable_moco.value() || !this->do_moco_)
        {
            GDEBUG_STREAM("Disable moco ...");
            moco = recon_res_.data_;
        }
        m1_sasha_moco->getObjectPtr()->data_ = moco;

        for (slc=0; slc<SLC; slc++)
        {
            for (s=0; s<S; s++)
            {
                for (n=0; n<N; n++)
                {
                    size_t ind = n + s*N + slc*N*S;
                    // ------------------------------
                    // Fill in headers
                    // ------------------------------
                    if (send_moco.value())
                    {
                        m1_sasha_moco->getObjectPtr()->headers_(n, s, slc) = recon_res_.headers_(n, s, slc);
                        m1_sasha_moco->getObjectPtr()->headers_(n, s, slc).image_index = (n + 1) + s*N + slc*N*S;
                        m1_sasha_moco->getObjectPtr()->headers_(n, s, slc).image_series_index += 100;
                        m1_sasha_moco->getObjectPtr()->headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                        m1_sasha_moco->getObjectPtr()->meta_[ind] = recon_res_.meta_[ind];
                        m1_sasha_moco->getObjectPtr()->meta_[ind].set(GADGETRON_DATA_ROLE, "GT");
                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_REGULAR);
                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_DATA_ROLE, "SASHA");
                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_DATA_ROLE, "MOCO");

                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_IMAGECOMMENT, "SASHA");
                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_IMAGECOMMENT, "MOCO");

                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_IMAGEPROCESSINGHISTORY, "MOCO");

                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA");
                        m1_sasha_moco->getObjectPtr()->meta_[ind].append(GADGETRON_SEQUENCEDESCRIPTION, "MOCO");

                        m1_sasha_moco->getObjectPtr()->meta_[ind].set(GADGETRON_IMAGENUMBER, (long)((n + 1) + s*N + slc*N*S));
                    }
                }
            }
        } // slc loop

        // ------------------------------
        // Enqueue
        // ------------------------------
        if (send_moco.value())
        {
            if (this->next()->putq(m1_sasha_moco) == -1)
            {
                GERROR("MoCoSashaGadget::process, passing m1_sasha_moco on to next gadget");
                return GADGET_FAIL;
            }
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    int MoCoSashaGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "MoCoSashaGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(MoCoSashaGadget)

}
