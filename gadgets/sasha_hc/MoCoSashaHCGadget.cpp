
#include "MoCoSashaHCGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "cmr_motion_correction.h"

namespace Gadgetron {

    MoCoSashaHCGadget::MoCoSashaHCGadget() : BaseClass()
    {
        this->do_moco_ = false;
    }

    MoCoSashaHCGadget::~MoCoSashaHCGadget()
    {
    }

    int MoCoSashaHCGadget::process_config(ACE_Message_Block* mb)
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
                    this->do_moco_ = true;
                    has_moco_in_proto = true;
                    GDEBUG_STREAM("MoCoSashaHCGadget, found MotionCorrection in protocol : " << this->do_moco_);
                }
            }
        }

        if (!has_moco_in_proto)
        {
            this->do_moco_ = true;
            GDEBUG_STREAM("MoCoSashaHCGadget, cannot find MotionCorrection in protocol : " << this->do_moco_);
        }

        return GADGET_OK;
    }

    int MoCoSashaHCGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("MoCoSashaHCGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "MoCoSashaHCGadget::process(...) starts ... ");

        process_called_times_++;

        std::stringstream os;
        os << "called_" << process_called_times_;
        std::string str = os.str();

        IsmrmrdImageArray* recon_res_ = m1->getObjectPtr();

        size_t encoding = (size_t)recon_res_->meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding < num_encoding_spaces_, GADGET_FAIL);

        std::string imageInfo;

        size_t RO = recon_res_->data_.get_size(0);
        size_t E1 = recon_res_->data_.get_size(1);
        size_t E2 = recon_res_->data_.get_size(2);
        size_t CHA = recon_res_->data_.get_size(3);
        size_t N = recon_res_->data_.get_size(4);
        size_t S = recon_res_->data_.get_size(5);
        size_t SLC = recon_res_->data_.get_size(6);

        GADGET_CHECK_RETURN(S == 2, GADGET_FAIL);

        size_t n, s, slc;

        // ----------------------------------------------------
        // compute diff images
        // ----------------------------------------------------
        hoNDArray<T> diff;
        diff.create(RO, E1, E2, CHA, N, 1, SLC);

        for (slc=0; slc<SLC; slc++)
        {
            T* pFirstN = &(recon_res_->data_(0, 0, 0, 0, 0, 0, slc));
            T* pSecondN = &(recon_res_->data_(0, 0, 0, 0, 0, 1, slc));

            hoNDArray<T> firstArray(RO, E1, E2, CHA, N, pFirstN);
            hoNDArray<T> secondArray(RO, E1, E2, CHA, N, pSecondN);

            T* pDiffN = &(diff(0, 0, 0, 0, 0, 0, slc));
            hoNDArray<T> diffArray(RO, E1, E2, CHA, N, pDiffN);

            // Subtract the magnitude images to avoid phase issues
        	hoNDArray<T> firstMag, secondMag;//, diff_mag, diff_mag_norm;
        	Gadgetron::abs(firstArray, firstMag);
        	Gadgetron::abs(secondArray, secondMag);
        	Gadgetron::subtract(firstMag, secondMag, diffArray);
        
        	Gadgetron::normalize(&diffArray);
        	diffArray *= (T)2048;
        	diffArray += (T)2048;
        }

        //for (slc = 0; slc < SLC; slc++)
        //{
        //    memcpy(&diff(0, 0, 0, 0, 0, 0, slc), &(recon_res_->data_(0, 0, 0, 0, 0, 1, slc)), sizeof(T) * RO * E1 * E2 * CHA * N);
        //}

        if (!debug_folder_full_path_.empty())
        {
            hoNDArray<T> data = recon_res_->data_;
            data.squeeze();
            gt_exporter_.export_array_complex(data, debug_folder_full_path_ + "MoCoSashaHC_data_" + str);
        }

        if (!debug_folder_full_path_.empty())
        {
            hoNDArray<T> data = diff;
            data.squeeze();
            gt_exporter_.export_array_complex(data, debug_folder_full_path_ + "MoCoSashaHC_diff_" + str);
        }

        // find the key frame and run the moco
        hoNDArray<real_value_type> diffMag;
        Gadgetron::abs(diff, diffMag);

        if (!debug_folder_full_path_.empty())
        {
            hoNDArray<real_value_type> data = diffMag;
            data.squeeze();
            gt_exporter_.export_array(data, debug_folder_full_path_ + "MoCoSashaHC_diffMag_" + str);
        }

        bool warp_input = false;

        hoNDArray<T> moco, diff_moco;
        moco.create(     RO, E1, E2, CHA, N, 1, SLC);
        diff_moco.create(RO, E1, E2, CHA, N, 1, SLC);

        for (slc=0; slc<SLC; slc++)
        {
            std::stringstream os;
            os << "slice_" << slc;
            std::string str_slc = os.str();

            real_value_type* pMag = &diffMag(0, 0, 0, 0, 0, 0, slc);

            hoNDArray<real_value_type> target(RO, E1, N, pMag);
            size_t key_frame(0);
            Gadgetron::find_key_frame_SSD_2DT(target, key_frame);

            // perform moco (on diff images)
            std::string moco_str = "GenericReconCartesianSpiritGadget, perform moco for " + str_slc;
            if (perform_timing.value()) { gt_timer_local_.start(moco_str.c_str()); }
            Gadgetron::hoImageRegContainer2DRegistration<ImageType, ImageType, double> reg;

            Gadgetron::perform_moco_fixed_key_frame_2DT(target, key_frame, (real_value_type)(10.0), iters_, bidirectional_moco.value(), warp_input, reg);

            if (perform_timing.value()) { gt_timer_local_.stop(); }

            // apply the deformation (on all images)
            moco_str = "GenericReconCartesianSpiritGadget, apply deformation field for " + str_slc;
            if (perform_timing.value()) { gt_timer_local_.start(moco_str.c_str()); }

            T* pFirstN = &(recon_res_->data_(0, 0, 0, 0, 0, 0, slc));
            hoNDArray<T> firstArray(RO, E1, E2, CHA, N, pFirstN);

            hoNDArray<real_value_type> firstArray_real, firstArray_imag;
            Gadgetron::complex_to_real_imag(firstArray, firstArray_real, firstArray_imag);

            firstArray_real.squeeze();
            firstArray_imag.squeeze();

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "called_" << process_called_times_ << "_slc" << slc;
                std::string str = os.str();

                gt_exporter_.export_array(firstArray_real, debug_folder_full_path_ + "MoCoSashaHC_firstArray_real_" + str);
                gt_exporter_.export_array(firstArray_imag, debug_folder_full_path_ + "MoCoSashaHC_firstArray_imag_" + str);
                gt_exporter_.export_array(target, debug_folder_full_path_ + "MoCoSashaHC_target_" + str);
            }

            hoNDArray<T> secondArray;
            secondArray.create(RO, E1, N);
            memcpy(secondArray.begin(), &(recon_res_->data_(0, 0, 0, 0, 0, 1, slc)), secondArray.get_number_of_bytes());
            hoNDArray<real_value_type> secondArray_real, secondArray_imag;
            Gadgetron::complex_to_real_imag(secondArray, secondArray_real, secondArray_imag);

            hoNDArray<double> dx, dy;
            reg.deformation_field_[0].to_NDArray(0, dx);
            reg.deformation_field_[1].to_NDArray(0, dy);

            hoNDArray<real_value_type> firstArray_real_moco, firstArray_imag_moco, target_moco;
            Gadgetron::apply_deformation_field(firstArray_real, dx, dy, firstArray_real_moco);
            Gadgetron::apply_deformation_field(firstArray_imag, dx, dy, firstArray_imag_moco);
            Gadgetron::apply_deformation_field(target,          dx, dy, target_moco);

            hoNDArray<real_value_type> secondArray_real_moco, secondArray_imag_moco;
            Gadgetron::apply_deformation_field(secondArray_real, dx, dy, secondArray_real_moco);
            Gadgetron::apply_deformation_field(secondArray_imag, dx, dy, secondArray_imag_moco);

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "called_" << process_called_times_ << "_slc" << slc;
                std::string str = os.str();

                gt_exporter_.export_array(firstArray_real_moco, debug_folder_full_path_ + "MoCoSashaHC_firstArray_real_moco_" + str);
                gt_exporter_.export_array(firstArray_imag_moco, debug_folder_full_path_ + "MoCoSashaHC_firstArray_imag_moco_" + str);
                gt_exporter_.export_array(target_moco, debug_folder_full_path_ + "MoCoSashaHC_target_moco_" + str);

                gt_exporter_.export_array(secondArray_real_moco, debug_folder_full_path_ + "MoCoSashaHC_secondArray_real_moco_" + str);
                gt_exporter_.export_array(secondArray_imag_moco, debug_folder_full_path_ + "MoCoSashaHC_secondArray_imag_moco_" + str);
            }

            hoNDArray<T> firstArrayMOCO;
            firstArrayMOCO.create(RO, E1, N);
            Gadgetron::real_imag_to_complex(firstArray_real_moco, firstArray_imag_moco, firstArrayMOCO);

            hoNDArray<T> secondArrayMOCO;
            secondArrayMOCO.create(RO, E1, N);
            Gadgetron::real_imag_to_complex(secondArray_real_moco, secondArray_imag_moco, secondArrayMOCO);

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "called_" << process_called_times_ << "_slc" << slc;
                std::string str = os.str();

                gt_exporter_.export_array_complex(firstArrayMOCO, debug_folder_full_path_ + "MoCoSashaHC_firstArrayMOCO_" + str);
                gt_exporter_.export_array_complex(secondArrayMOCO, debug_folder_full_path_ + "MoCoSashaHC_secondArrayMOCO_" + str);
            }

            {
                hoNDArray<real_value_type> target2;
                Gadgetron::abs(secondArrayMOCO, target2);
                Gadgetron::hoImageRegContainer2DRegistration<ImageType, ImageType, double> reg2;
                Gadgetron::perform_moco_fixed_key_frame_2DT(target2, key_frame, (real_value_type)(regularization_hilbert_strength.value()), iters_, bidirectional_moco.value(), warp_input, reg2);

                hoNDArray<double> dx2, dy2;
                reg2.deformation_field_[0].to_NDArray(0, dx2);
                reg2.deformation_field_[1].to_NDArray(0, dy2);

                hoNDArray<real_value_type> firstArray_real2, firstArray_imag2;
                Gadgetron::complex_to_real_imag(firstArrayMOCO, firstArray_real2, firstArray_imag2);

                hoNDArray<real_value_type> firstArray_real_moco2, firstArray_imag_moco2, target_moco2;
                Gadgetron::apply_deformation_field(firstArray_real2, dx2, dy2, firstArray_real_moco2);
                Gadgetron::apply_deformation_field(firstArray_imag2, dx2, dy2, firstArray_imag_moco2);
                Gadgetron::apply_deformation_field(target2, dx2, dy2, target_moco2);

                hoNDArray<T> firstArrayMOCO2;
                firstArrayMOCO2.create(RO, E1, N);
                Gadgetron::real_imag_to_complex(firstArray_real_moco2, firstArray_imag_moco2, firstArrayMOCO2);

                /*std::vector<double> mean_deform, max_deform, mean_log_jac, max_log_jac;
                Gadgetron::compute_deformation_jacobian(dx2, dy2, mean_deform, max_deform, mean_log_jac, max_log_jac);*/

                memcpy(&moco(0, 0, 0, 0, 0, 0, slc), firstArrayMOCO2.begin(), firstArrayMOCO2.get_number_of_bytes());

                if (!debug_folder_full_path_.empty())
                {
                    std::stringstream os;
                    os << "called_" << process_called_times_ << "_slc" << slc;
                    std::string str = os.str();

                    gt_exporter_.export_array(target2, debug_folder_full_path_ + "MoCoSashaHC_target2_" + str);
                    gt_exporter_.export_array(target_moco2, debug_folder_full_path_ + "MoCoSashaHC_target_moco2_" + str);

                    gt_exporter_.export_array(dx2, debug_folder_full_path_ + "MoCoSashaHC_dx2_" + str);
                    gt_exporter_.export_array(dy2, debug_folder_full_path_ + "MoCoSashaHC_dy2_" + str);
                    gt_exporter_.export_array_complex(firstArrayMOCO2, debug_folder_full_path_ + "MoCoSashaHC_firstArrayMOCO2_" + str);
                }
            }

            pFirstN = &(diff_moco(0, 0, 0, 0, 0, 0, slc));
            hoNDArray<T> diffMOCO(RO, E1, N, pFirstN);
            diffMOCO.copyFrom(target_moco);

            if (perform_timing.value()) { gt_timer_local_.stop(); }

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "called_" << process_called_times_ << "_slc" << slc;
                std::string str = os.str();

                hoNDArray<T> data;
                data = diffMOCO; data.squeeze();
                gt_exporter_.export_array_complex(data,       debug_folder_full_path_ + "MoCoSashaHC_diffMag_MOCO_"    + str);

                gt_exporter_.export_array(dx, debug_folder_full_path_ + "MoCoSashaHC_dx_" + str);
                gt_exporter_.export_array(dy, debug_folder_full_path_ + "MoCoSashaHC_dy_" + str);
            }
        }

        if (!debug_folder_full_path_.empty())
        {
            std::stringstream os;
            os << "called_" << process_called_times_ << "_slc" << slc;
            std::string str = os.str();

            hoNDArray<T> data;
            data = moco; data.squeeze();
            gt_exporter_.export_array_complex(data, debug_folder_full_path_ + "MoCoSashaHC_moco_res_" + str);
        }

        GDEBUG_CONDITION_STREAM(verbose.value(), "MoCoSashaHCGadget::process(...) ends ... ");

		// ----------------------------------------------------
		// format data to be passed on
		// ----------------------------------------------------
		Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha           = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
		Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_moco      = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
		Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_hc        = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
		Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_diff      = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
		Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_diff_moco = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();

		if (send_ori.value())
		{
			m1_sasha->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
			m1_sasha->getObjectPtr()->headers_.create(N, 1, SLC);
			m1_sasha->getObjectPtr()->meta_.resize(N*1*SLC);
		}

		if (send_moco.value())
		{
			m1_sasha_moco->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
			m1_sasha_moco->getObjectPtr()->headers_.create(N, 1, SLC);
			m1_sasha_moco->getObjectPtr()->meta_.resize(N*1*SLC);
		}

		if (send_hc.value())
		{
			m1_sasha_hc->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
			m1_sasha_hc->getObjectPtr()->headers_.create(N, 1, SLC);
			m1_sasha_hc->getObjectPtr()->meta_.resize(N*1*SLC);
		}

		if (send_diff.value())
		{
			m1_sasha_diff->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
			m1_sasha_diff->getObjectPtr()->headers_.create(N, 1, SLC);
			m1_sasha_diff->getObjectPtr()->meta_.resize(N*1*SLC);
		}

		if (send_diff.value() && send_moco.value())
		{
			m1_sasha_diff_moco->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
			m1_sasha_diff_moco->getObjectPtr()->headers_.create(N, 1, SLC);
			m1_sasha_diff_moco->getObjectPtr()->meta_.resize(N*1*SLC);
		}

        if(this->disable_moco.value() || !this->do_moco_)
        {
            GDEBUG_STREAM("Disable moco ...");
            moco = recon_res_->data_;
        }

		for (slc=0; slc<SLC; slc++)
		{
			for (n=0; n<N; n++)
			{
				// ------------------------------
				// Copy image data
				// ------------------------------
				if (send_ori.value())
				{
					memcpy( &(m1_sasha->getObjectPtr()->data_(     0, 0, 0, 0, n, 0, slc)), 
					        &(recon_res_->data_(                   0, 0, 0, 0, n, 0, slc)), sizeof(T)*RO*E1*E2*CHA);
				}

				if (send_moco.value())
				{
					memcpy( &(m1_sasha_moco->getObjectPtr()->data_(0, 0, 0, 0, n, 0, slc)),
					        &(moco(                                0, 0, 0, 0, n, 0, slc)), sizeof(T)*RO*E1*E2*CHA);
				}

				if (send_hc.value())
				{
					memcpy( &(m1_sasha_hc->getObjectPtr()->data_(  0, 0, 0, 0, n, 0, slc)), 
					        &(recon_res_->data_(                   0, 0, 0, 0, n, 1, slc)), sizeof(T)*RO*E1*E2*CHA);
				}

				if (send_diff.value())
				{
					memcpy( &(m1_sasha_diff->getObjectPtr()->data_(0, 0, 0, 0, n, 0, slc)), 
					        &(diff(                                0, 0, 0, 0, n, 0, slc)), sizeof(T)*RO*E1*E2*CHA);
				}

				if (send_diff.value() && send_moco.value())
				{
					memcpy( &(m1_sasha_diff_moco->getObjectPtr()->data_(0, 0, 0, 0, n, 0, slc)), 
					        &(diff_moco(                                0, 0, 0, 0, n, 0, slc)), sizeof(T)*RO*E1*E2*CHA);						
				}

				// TS time
				std::stringstream ss_ts_time;
				ss_ts_time << "TS" << (int16_t)recon_res_->headers_(n, 0, slc).user_int[4];
//				GDEBUG_STREAM("n: " << n << " " << ss_ts_time.str());

				// ------------------------------
				// Fill in headers
				// ------------------------------
				if (send_ori.value())
				{
                    m1_sasha->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 0, slc);
                    m1_sasha->getObjectPtr()->headers_(n, 0, slc).image_index = (n + 1) + N * slc;
                    m1_sasha->getObjectPtr()->headers_(n, 0, slc).image_series_index += 0;
					m1_sasha->getObjectPtr()->headers_(n, 0, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                    m1_sasha->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+0*N+slc*N*S];
                    m1_sasha->getObjectPtr()->meta_[n + slc * N].set(GADGETRON_DATA_ROLE, "GT");
                    m1_sasha->getObjectPtr()->meta_[n + slc * N].append(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_REGULAR);
                    m1_sasha->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA");

                    m1_sasha->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA");
                    m1_sasha->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA");

					m1_sasha->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, ss_ts_time.str().c_str());
                    m1_sasha->getObjectPtr()->meta_[n+slc*N].append("Skip_processing_after_recon", "true");

                    m1_sasha->getObjectPtr()->meta_[n + slc * N].set(GADGETRON_IMAGENUMBER, (long)((n + 1) + N * slc));
				}

				if (send_moco.value())
				{
                    m1_sasha_moco->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 0, slc);
                    m1_sasha_moco->getObjectPtr()->headers_(n, 0, slc).image_index = (n + 1) + N * slc;
                    m1_sasha_moco->getObjectPtr()->headers_(n, 0, slc).image_series_index += 100;
					m1_sasha_moco->getObjectPtr()->headers_(n, 0, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+0*N+slc*N*S];
                    m1_sasha_moco->getObjectPtr()->meta_[n + slc * N].set(GADGETRON_DATA_ROLE, "GT");
                    m1_sasha_moco->getObjectPtr()->meta_[n + slc * N].append(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_REGULAR);
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA");
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "MOCO");

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA");
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "MOCO");

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGEPROCESSINGHISTORY, "MOCO");

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA");
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "MOCO");

                    m1_sasha_moco->getObjectPtr()->meta_[n + slc * N].set(GADGETRON_IMAGENUMBER, (long)((n + 1) + N * slc));
				}

				if (send_hc.value())
				{
                    m1_sasha_hc->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 1, slc);
                    m1_sasha_hc->getObjectPtr()->headers_(n, 0, slc).image_series_index += 1; // Not +2 because the second contrast already gave it a +1
					m1_sasha_hc->getObjectPtr()->headers_(n, 0, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+1*N+slc*N*S];
                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA_HC");

                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA_HC");
                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA_HC");
                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N].append("Skip_processing_after_recon", "true");
				}

				if (send_diff.value())
				{
                    m1_sasha_diff->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 0, slc);
                    m1_sasha_diff->getObjectPtr()->headers_(n, 0, slc).image_series_index += 3;
					m1_sasha_diff->getObjectPtr()->headers_(n, 0, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+0*N+slc*N*S];
                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE,           "SASHA_HC_DIFF");
                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT,        "SASHA_HC_DIFF");
                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA_HC_DIFF");
                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append("Skip_processing_after_recon", "true");
                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append("Skip_processing_after_recon", "true");
				}

				if (send_diff.value() && send_moco.value())
				{
                    m1_sasha_diff_moco->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 0, slc);
                    m1_sasha_diff_moco->getObjectPtr()->headers_(n, 0, slc).image_series_index += 4;
					m1_sasha_diff_moco->getObjectPtr()->headers_(n, 0, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+0*N+slc*N*S];
                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA_HC_DIFF");
                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "MOCO");

                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA_HC_DIFF");
                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "MOCO");

                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGEPROCESSINGHISTORY, "MOCO");

                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA_HC_DIFF");
                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "MOCO");
                    m1_sasha_diff_moco->getObjectPtr()->meta_[n+slc*N].append("Skip_processing_after_recon", "true");
				}
			}
		} // slc loop

        // ------------------------------
        // Enqueue
        // ------------------------------
        if (send_ori.value())
        {
            if (this->next()->putq(m1_sasha) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha on to next gadget");
                return GADGET_FAIL;
            }
        }

        if (send_moco.value())
        {
            if (this->next()->putq(m1_sasha_moco) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha_moco on to next gadget");
                return GADGET_FAIL;
            }
        }

        if (send_hc.value())
        {
            if (this->next()->putq(m1_sasha_hc) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha_hc on to next gadget");
                return GADGET_FAIL;
            }
        }

        if (send_diff.value())
        {
            if (this->next()->putq(m1_sasha_diff) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha_diff on to next gadget");
                return GADGET_FAIL;
            }
        }

        if (send_diff.value() && send_moco.value())
        {
            if (this->next()->putq(m1_sasha_diff_moco) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha_diff_moco on to next gadget");
                return GADGET_FAIL;
            }
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    int MoCoSashaHCGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "MoCoSashaHCGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(MoCoSashaHCGadget)

}
