
#include "GtImageMoCoAvePSIRGadget.h"

namespace Gadgetron {

    GtImageMoCoAvePSIRGadget::GtImageMoCoAvePSIRGadget() : BaseClass()
    {
        PD_set_ = 1;
        apply_PD_filtering_ = true;;
        filter_width_ = 5;
        scale_factor_after_SCC_ = 250;
        offset_after_SCC_ = 0;
        windowing_high_end_percentile_ = 0.95;

        scc_strategy_ = "FFDM";
        num_of_refinement_FFD_ = 2;
        num_of_refinement_max_FFD_ = 7;

        preserve_PD_for_scc_ = true;

        noise_masking_ = true;

        thres_ratio_noise_masking_ = 3.0f;

        perform_scc_PSIR_ = true;
        perform_scc_mag_IR_ = false;

        send_ori_PSIR_ = false;
        send_ori_mag_IR_ = false;
        send_ori_mag_PD_ = false;

        send_moco_PSIR_ = false;
        send_moco_mag_IR_ = false;
        send_moco_mag_PD_ = false;

        send_moco_ave_PSIR_ = true;
        send_moco_ave_mag_IR_ = true;
        send_moco_ave_mag_PD_ = false;

        send_no_scc_mag_IR_ = false;
        send_no_scc_PSIR_ = false;
    }

    GtImageMoCoAvePSIRGadget::~GtImageMoCoAvePSIRGadget()
    {
    }

    bool GtImageMoCoAvePSIRGadget::readParameters()
    {
        try
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "------> GtImageMoCoAvePSIRGadget parameters <------");

            debugFolder_PSIR_ = debugFolder_PSIR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "debugFolder_PSIR_ is " << debugFolder_PSIR_);

            if (!debugFolder_PSIR_.empty())
            {
                Gadgetron::get_debug_folder_path(debugFolder_PSIR_, debugFolder_PSIR_fullPath_);
            }
            else
            {
                GDEBUG_STREAM("GtImageMoCoAvePSIRGadget, debugFolder_PSIR is not set ...");
            }

            PD_set_ = PD_set.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "PD_set_ is " << PD_set_);

            apply_PD_filtering_ = apply_PD_filtering.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "apply_PD_filtering_ is " << apply_PD_filtering_);

            preserve_PD_for_scc_ = preserve_PD_for_scc.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "preserve_PD_for_scc_ is " << preserve_PD_for_scc_);

            filter_width_ = filter_width.value();
            if (filter_width_ <= 0) filter_width_ = 5;
            GDEBUG_CONDITION_STREAM(verbose.value(), "filter_width_ is " << filter_width_);

            scale_factor_after_SCC_ = scale_factor_after_SCC.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "scale_factor_after_SCC_ is " << scale_factor_after_SCC_);

            offset_after_SCC_ = offset_after_SCC.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "offset_after_SCC_ is " << offset_after_SCC_);

            windowing_high_end_percentile_ = windowing_high_end_percentile.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "windowing_high_end_percentile_ is " << windowing_high_end_percentile_);

            scc_strategy_ = scc_strategy.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "scc_strategy_ is " << scc_strategy_);

            num_of_refinement_FFD_ = num_of_refinement_FFD.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "num_of_refinement_FFD_ is " << num_of_refinement_FFD_);

            num_of_refinement_max_FFD_ = num_of_refinement_max_FFD.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "num_of_refinement_max_FFD_ is " << num_of_refinement_max_FFD_);

            noise_masking_ = noise_masking.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "noise_masking_ is " << noise_masking_);

            thres_ratio_noise_masking_ = thres_ratio_noise_masking.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "thres_ratio_noise_masking_ is " << thres_ratio_noise_masking_);

            perform_scc_PSIR_ = perform_scc_PSIR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "perform_scc_PSIR_ is " << perform_scc_PSIR_);

            perform_scc_mag_IR_ = perform_scc_mag_IR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "perform_scc_mag_IR_ is " << perform_scc_mag_IR_);

            // ------------------------------------------------------------

            send_ori_PSIR_ = send_ori_PSIR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_ori_PSIR_ is " << send_ori_PSIR_);

            send_ori_mag_IR_ = send_ori_mag_IR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_ori_mag_IR_ is " << send_ori_mag_IR_);

            send_ori_mag_PD_ = send_ori_mag_PD.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_ori_mag_PD_ is " << send_ori_mag_PD_);

            // ------------------------------------------------------------

            send_moco_PSIR_ = send_moco_PSIR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_moco_PSIR_ is " << send_moco_PSIR_);

            send_moco_mag_IR_ = send_moco_mag_IR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_moco_mag_IR_ is " << send_moco_mag_IR_);

            send_moco_mag_PD_ = send_moco_mag_PD.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_moco_mag_PD_ is " << send_moco_mag_PD_);

            // ------------------------------------------------------------

            send_moco_ave_PSIR_ = send_moco_ave_PSIR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_moco_ave_PSIR_ is " << send_moco_ave_PSIR_);

            send_moco_ave_mag_IR_ = send_moco_ave_mag_IR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_moco_ave_mag_IR_ is " << send_moco_ave_mag_IR_);

            send_moco_ave_mag_PD_ = send_moco_ave_mag_PD.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_moco_ave_mag_PD_ is " << send_moco_ave_mag_PD_);

            GDEBUG_CONDITION_STREAM(verbose.value(), "-----------------------------------------------");

            send_no_scc_mag_IR_ = send_no_scc_mag_IR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_no_scc_mag_IR_ is " << send_no_scc_mag_IR_);

            send_no_scc_PSIR_ = send_no_scc_PSIR.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_no_scc_PSIR_ is " << send_no_scc_PSIR_);
        }
        catch (...)
        {
            GERROR_STREAM("Errors in GtImageMoCoAvePSIRGadget::readParameters() ... ");
            return false;
        }

        return true;
    }

    int GtImageMoCoAvePSIRGadget::process_config(ACE_Message_Block* mb)
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
            throw;
            return GADGET_FAIL;
        }

        if (h.sequenceParameters.is_present())
        {
            if (h.sequenceParameters.get().TI.is_present())
            {
                TI_ = h.sequenceParameters.get().TI.get();
            }
        }
        else
        {
            GWARN_STREAM("Inversion time does not exist in the seq protocols ... ");
        }

        GADGET_CHECK_RETURN(this->readParameters(), GADGET_FAIL);

        psir_reconer_.debugFolder_ = this->debugFolder_PSIR_fullPath_;
        psir_reconer_.verbose_ = this->verbose.value();
        psir_reconer_.performTiming_ = this->perform_timing.value();

        psir_reconer_.row_PD_ = PD_set_;
        psir_reconer_.col_compute_PSIR_windowing_ = 0;
        psir_reconer_.perform_SCC_PSIR_ = this->perform_scc_PSIR_;
        psir_reconer_.perform_SCC_MagIR_ = this->perform_scc_mag_IR_;
        psir_reconer_.apply_PD_filtering_ = this->apply_PD_filtering_;
        psir_reconer_.preserve_PD_for_scc_ = this->preserve_PD_for_scc_;
        psir_reconer_.scale_factor_after_SCC_ = this->scale_factor_after_SCC_;
        psir_reconer_.offset_after_SCC_ = this->offset_after_SCC_;
        psir_reconer_.windowing_high_end_percentile_ = this->windowing_high_end_percentile_;

        psir_reconer_.filter_width_[0] = this->filter_width_;
        psir_reconer_.filter_width_[1] = this->filter_width_;

        if (inteprolation_on_)
        {
            GDEBUG_STREAM("Interpolation is ON; double the surface coil correction filter width ... ");
            psir_reconer_.filter_width_[0] *= 2;
            psir_reconer_.filter_width_[1] *= 2;
        }

        psir_reconer_.scc_strategy_ = this->scc_strategy_;
        psir_reconer_.num_of_refinement_FFD_ = this->num_of_refinement_FFD_;
        psir_reconer_.num_of_refinement_max_FFD_ = this->num_of_refinement_max_FFD_;
        psir_reconer_.noise_masking_ = this->noise_masking_;
        psir_reconer_.thres_ratio_noise_masking_ = this->thres_ratio_noise_masking_;

        psir_reconer_.perform_PSIR_ = true;

        return GADGET_OK;
    }

    int GtImageMoCoAvePSIRGadget::processImageBuffer(ImageBufferType& ori)
    {
        GDEBUG_CONDITION_STREAM(verbose.value(), "GtImageMoCoAvePSIRGadget::processImageBuffer(...) starts ... ");

        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        std::vector<size_t> dims = ori.get_dimensions();
        GDEBUG_CONDITION_STREAM(verbose.value(), "[Cha Slice Con Phase Rep Set Ave] = [" << dims[0] << " " << dims[1] << " " << dims[2] << " " << dims[3] << " " << dims[4] << " " << dims[5] << " " << dims[6] << "]");

        if (ori(0).get_number_of_elements() > 0)
        {
            try
            {
                psir_reconer_.intensity_scale_factor_ = (float)ori(0).attrib_.as_double(GADGETRON_IMAGE_SCALE_RATIO, 0);
                GDEBUG_CONDITION_STREAM(verbose.value(), "image intensity scaling factor is found to be " << psir_reconer_.intensity_scale_factor_);
            }
            catch (...)
            {
                GWARN_STREAM("image intensity scaling factor can not be found ... ");
                psir_reconer_.intensity_scale_factor_ = 1;
            }
        }

        // --------------------------------------------------------------------------------
        // ori
        // --------------------------------------------------------------------------------
        if (send_ori_)
        {
            GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_, processStr, dataRole), GADGET_FAIL);
        }

        if (send_ori_mag_IR_ || send_ori_PSIR_ || send_ori_mag_PD_)
        {
            ImageBufferType oriMagIR, oriMagIRNoSCC, oriPSIR, oriPSIRNoSCC, oriMagPD;
            std::vector<float> oriWindowCenter, oriWindowWidth;

            GADGET_CHECK_RETURN(this->performPSIR(ori, oriMagIR, oriMagIRNoSCC, oriPSIR, oriPSIRNoSCC, oriMagPD, oriWindowCenter, oriWindowWidth), GADGET_FAIL);

            processStr.clear();

            if (send_ori_mag_IR_)
            {
                if (this->perform_scc_mag_IR_) processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);

                dataRole.clear();
                if (this->perform_scc_mag_IR_) dataRole.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                dataRole.push_back(GADGETRON_IMAGE_MAGIR);

                GADGET_CHECK_RETURN(this->sendOutImages(oriMagIR, image_series_num_ + 1, processStr, dataRole), GADGET_FAIL);

                if (this->perform_scc_mag_IR_ && this->send_no_scc_mag_IR_)
                {
                    processStr.clear();
                    dataRole.clear();

                    dataRole.push_back(GADGETRON_IMAGE_MAGIR);
                    GADGET_CHECK_RETURN(this->sendOutImages(oriMagIRNoSCC, image_series_num_ + 12, processStr, dataRole), GADGET_FAIL);
                }

                if (this->perform_scc_PSIR_ && this->send_no_scc_mag_IR_)
                {
                    processStr.clear();
                    dataRole.clear();

                    dataRole.push_back(GADGETRON_IMAGE_PSIR);
                    GADGET_CHECK_RETURN(this->sendOutImages(oriPSIRNoSCC, image_series_num_ + 15, processStr, dataRole), GADGET_FAIL);
                }
            }

            if (send_ori_mag_PD_)
            {
                dataRole.clear();
                dataRole.resize(2);
                dataRole[0] = GADGETRON_IMAGE_PD;

                processStr.clear();

                GADGET_CHECK_RETURN(this->sendOutImages(oriMagPD, image_series_num_ + 2, processStr, dataRole), GADGET_FAIL);
            }

            if (send_ori_PSIR_)
            {
                processStr.clear();
                if (this->perform_scc_PSIR_) processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                processStr.push_back(GADGETRON_IMAGE_PSIR);

                dataRole.clear();
                if (this->perform_scc_PSIR_) dataRole.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                dataRole.push_back(GADGETRON_IMAGE_PSIR);

                GADGET_CHECK_RETURN(this->sendOutImages(oriPSIR, image_series_num_ + 3, processStr, dataRole, oriWindowCenter, oriWindowWidth), GADGET_FAIL);
            }

            GADGET_CHECK_RETURN(this->releaseImageBuffer(oriMagIR), GADGET_FAIL);
            GADGET_CHECK_RETURN(this->releaseImageBuffer(oriMagIRNoSCC), GADGET_FAIL);
            GADGET_CHECK_RETURN(this->releaseImageBuffer(oriPSIRNoSCC), GADGET_FAIL);
            GADGET_CHECK_RETURN(this->releaseImageBuffer(oriPSIR), GADGET_FAIL);
            GADGET_CHECK_RETURN(this->releaseImageBuffer(oriMagPD), GADGET_FAIL);
        }

        ImageBufferType moco, ave;
        if (moco_ave_)
        {
            bool mocoPerformed(true);
            GADGET_CHECK_RETURN(this->performMoCoAveraging(ori, moco, ave, mocoPerformed), GADGET_FAIL);

            // --------------------------------------------------------------------------------
            // moco
            // --------------------------------------------------------------------------------
            if (send_moco_)
            {
                processStr.clear();
                if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);

                dataRole.clear();
                if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                GADGET_CHECK_RETURN(this->sendOutImages(moco, image_series_num_ + 4, processStr, dataRole), GADGET_FAIL);
            }

            if (send_moco_mag_IR_ || send_moco_PSIR_ || send_moco_mag_PD_)
            {
                ImageBufferType mocoMagIR, mocoMagIRNoSCC, mocoPSIR, mocoPSIRNoSCC, mocoMagPD;
                std::vector<float> mocoWindowCenter, mocoWindowWidth;

                GADGET_CHECK_RETURN(this->performPSIR(moco, mocoMagIR, mocoMagIRNoSCC, mocoPSIR, mocoPSIRNoSCC, mocoMagPD, mocoWindowCenter, mocoWindowWidth), GADGET_FAIL);

                if (send_moco_mag_IR_)
                {
                    processStr.clear();
                    if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                    if (this->perform_scc_mag_IR_) processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);

                    dataRole.clear();
                    if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                    if (this->perform_scc_mag_IR_) dataRole.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                    dataRole.push_back(GADGETRON_IMAGE_MAGIR);

                    GADGET_CHECK_RETURN(this->sendOutImages(mocoMagIR, image_series_num_ + 5, processStr, dataRole), GADGET_FAIL);
                    processStr.pop_back();

                    if (this->perform_scc_mag_IR_ && this->send_no_scc_mag_IR_)
                    {
                        processStr.clear();
                        dataRole.clear();

                        if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                        if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                        dataRole.push_back(GADGETRON_IMAGE_MAGIR);

                        GADGET_CHECK_RETURN(this->sendOutImages(mocoMagIRNoSCC, image_series_num_ + 13, processStr, dataRole), GADGET_FAIL);
                    }

                    if (this->perform_scc_mag_IR_ && this->send_no_scc_PSIR_)
                    {
                        processStr.clear();
                        dataRole.clear();

                        if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                        if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                        dataRole.push_back(GADGETRON_IMAGE_PSIR);

                        GADGET_CHECK_RETURN(this->sendOutImages(mocoPSIRNoSCC, image_series_num_ + 16, processStr, dataRole), GADGET_FAIL);
                    }
                }

                if (send_moco_mag_PD_)
                {
                    processStr.clear();
                    if (mocoPerformed)  processStr.push_back(GADGETRON_IMAGE_MOCO);

                    dataRole.clear();
                    if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                    dataRole.push_back(GADGETRON_IMAGE_PD);

                    GADGET_CHECK_RETURN(this->sendOutImages(mocoMagPD, image_series_num_ + 6, processStr, dataRole), GADGET_FAIL);
                }

                if (send_moco_PSIR_)
                {
                    processStr.clear();
                    if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                    if (this->perform_scc_PSIR_) processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                    processStr.push_back(GADGETRON_IMAGE_PSIR);

                    dataRole.clear();
                    if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                    if (this->perform_scc_PSIR_) dataRole.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                    dataRole.push_back(GADGETRON_IMAGE_PSIR);

                    GADGET_CHECK_RETURN(this->sendOutImages(mocoPSIR, image_series_num_ + 7, processStr, dataRole, mocoWindowCenter, mocoWindowWidth), GADGET_FAIL);
                    processStr.pop_back();
                }

                GADGET_CHECK_RETURN(this->releaseImageBuffer(mocoMagIR), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(mocoMagIRNoSCC), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(mocoPSIRNoSCC), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(mocoPSIR), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(mocoMagPD), GADGET_FAIL);
            }

            // --------------------------------------------------------------------------------
            // moco ave
            // --------------------------------------------------------------------------------
            if (send_moco_ave_)
            {
                processStr.clear();
                if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                processStr.push_back(GADGETRON_IMAGE_AVE);

                if (send_moco_ave_)
                {
                    dataRole.clear();
                    if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                    dataRole.push_back(GADGETRON_IMAGE_AVE);
                    GADGET_CHECK_RETURN(this->sendOutImages(ave, image_series_num_ + 8, processStr, dataRole), GADGET_FAIL);
                }
            }

            if (send_moco_ave_mag_IR_ || send_moco_ave_PSIR_ || send_moco_ave_mag_PD_)
            {
                ImageBufferType aveMagIR, aveMagIRNoSCC, avePSIR, avePSIRNoSCC, aveMagPD;
                std::vector<float> aveWindowCenter, aveWindowWidth;

                GADGET_CHECK_RETURN(this->performPSIR(ave, aveMagIR, aveMagIRNoSCC, avePSIR, avePSIRNoSCC, aveMagPD, aveWindowCenter, aveWindowWidth), GADGET_FAIL);

                if (send_moco_ave_mag_IR_)
                {
                    processStr.clear();
                    if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                    processStr.push_back(GADGETRON_IMAGE_AVE);
                    if (this->perform_scc_mag_IR_) processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);

                    dataRole.clear();
                    if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                    dataRole.push_back(GADGETRON_IMAGE_AVE);
                    if (this->perform_scc_mag_IR_) dataRole.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                    dataRole.push_back(GADGETRON_IMAGE_MAGIR);

                    GADGET_CHECK_RETURN(this->sendOutImages(aveMagIR, image_series_num_ + 9, processStr, dataRole), GADGET_FAIL);
                    processStr.pop_back();

                    if (this->perform_scc_mag_IR_ && this->send_no_scc_mag_IR_)
                    {
                        processStr.clear();
                        if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                        processStr.push_back(GADGETRON_IMAGE_AVE);

                        dataRole.clear();
                        if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                        dataRole.push_back(GADGETRON_IMAGE_AVE);
                        dataRole.push_back(GADGETRON_IMAGE_MAGIR);

                        GADGET_CHECK_RETURN(this->sendOutImages(aveMagIRNoSCC, image_series_num_ + 14, processStr, dataRole), GADGET_FAIL);
                    }

                    if (this->perform_scc_mag_IR_ && this->send_no_scc_PSIR_)
                    {
                        processStr.clear();
                        if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                        processStr.push_back(GADGETRON_IMAGE_AVE);

                        dataRole.clear();
                        if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                        dataRole.push_back(GADGETRON_IMAGE_AVE);
                        dataRole.push_back(GADGETRON_IMAGE_PSIR);

                        GADGET_CHECK_RETURN(this->sendOutImages(avePSIRNoSCC, image_series_num_ + 17, processStr, dataRole), GADGET_FAIL);
                    }
                }

                if (send_moco_ave_mag_PD_)
                {
                    processStr.clear();
                    if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                    processStr.push_back(GADGETRON_IMAGE_AVE);

                    dataRole.clear();
                    if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                    dataRole.push_back(GADGETRON_IMAGE_AVE);
                    dataRole.push_back(GADGETRON_IMAGE_PD);

                    GADGET_CHECK_RETURN(this->sendOutImages(aveMagPD, image_series_num_ + 10, processStr, dataRole), GADGET_FAIL);
                }

                if (send_moco_ave_PSIR_)
                {
                    processStr.clear();
                    if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                    processStr.push_back(GADGETRON_IMAGE_AVE);
                    if (this->perform_scc_PSIR_) processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                    processStr.push_back(GADGETRON_IMAGE_PSIR);

                    dataRole.clear();
                    if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                    dataRole.push_back(GADGETRON_IMAGE_AVE);
                    if (this->perform_scc_PSIR_) dataRole.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
                    dataRole.push_back(GADGETRON_IMAGE_PSIR);

                    GADGET_CHECK_RETURN(this->sendOutImages(avePSIR, image_series_num_ + 11, processStr, dataRole, aveWindowCenter, aveWindowWidth), GADGET_FAIL);
                    processStr.pop_back();
                }

                GADGET_CHECK_RETURN(this->releaseImageBuffer(aveMagIR), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(aveMagIRNoSCC), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(avePSIRNoSCC), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(avePSIR), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(aveMagPD), GADGET_FAIL);
            }
        }

        // release the moco images
        GADGET_CHECK_RETURN(this->releaseImageBuffer(moco), GADGET_FAIL);

        // release the moco ave images
        GADGET_CHECK_RETURN(this->releaseImageBuffer(ave), GADGET_FAIL);

        GDEBUG_CONDITION_STREAM(verbose.value(), "GtImageMoCoAvePSIRGadget::process(...) ends ... ");

        return GADGET_OK;
    }

    bool GtImageMoCoAvePSIRGadget::performPSIR(ImageBufferType& input, ImageBufferType& magIR, ImageBufferType& magIRNoSCC, ImageBufferType& psIR, ImageBufferType& psIRNoSCC, ImageBufferType& magPD, std::vector<float>& wCenter, std::vector<float>& wWidth)
    {
        size_t CHA = input.get_size(0);
        size_t SLC = input.get_size(1);
        size_t CON = input.get_size(2);
        size_t PHS = input.get_size(3);
        size_t REP = input.get_size(4);
        size_t SET = input.get_size(5);
        size_t AVE = input.get_size(6);

        std::vector<size_t> dim;
        input.get_dimensions(dim);

        dim[5] = 1; // only one set is outputted

        magIR.create(dim);
        magIRNoSCC.create(dim);
        psIR.create(dim);
        psIRNoSCC.create(dim);
        magPD.create(dim);

        size_t cha, slc, con, phs, rep, set, ave;

        ImageContainer2DType inputContainer, magIRContainer, magIRNoSCCContainer, PSIRContainer, psIRNoSCCContainer, magPDContainer;
        float windowCenter = -1;
        float windowWidth = -1;
        psir_reconer_.compute_PSIR_windowing_ = true;

        wCenter.resize(SLC, -1);
        wWidth.resize(SLC, -1);

        std::vector<size_t> cols(SET, CON);
        inputContainer.create(cols, false);

        for (cha = 0; cha < CHA; cha++)
        {
            for (ave = 0; ave < AVE; ave++)
            {
                for (slc = 0; slc < SLC; slc++)
                {
                    windowCenter = -1;
                    windowWidth = -1;

                    psir_reconer_.compute_PSIR_windowing_ = true;

                    ImageType gmap;
                    bool gMapFound = this->findGFactorMap(slc, gmap);
                    float gMap_scale_factor = 1.0f;
                    if (gMapFound)
                    {
                        GDEBUG_CONDITION_STREAM(verbose.value(), "Gfactor map for slc - " << slc << " is found ... ");

                        if (gmap.attrib_.length(GADGETRON_IMAGE_SCALE_RATIO) > 0)
                        {
                            gMap_scale_factor = (float)gmap.attrib_.as_double(GADGETRON_IMAGE_SCALE_RATIO, 0);
                            GDEBUG_CONDITION_STREAM(verbose.value(), "Gfactor map for slc - " << slc << " has scale factor " << gMap_scale_factor);

                            Gadgetron::scal((float)(1.0 / gMap_scale_factor), gmap);
                        }

                        if (!debug_folder_full_path_.empty())
                        {
                            std::ostringstream ostr;
                            ostr << "gfactor_map_slc" << slc;
                            if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(gmap, debug_folder_full_path_ + ostr.str());
                        }
                    }
                    else
                    {
                        GDEBUG_CONDITION_STREAM(true, "Gfactor map for slc - " << slc << " is NOT found ... ");
                    }

                    for (phs = 0; phs < PHS; phs++)
                    {
                        for (rep = 0; rep < REP; rep++)
                        {
                            // compute PSIR for every SET/CON

                            for (con = 0; con < CON; con++)
                            {
                                for (set = 0; set < SET; set++)
                                {
                                    inputContainer.set(&input(cha, slc, con, phs, rep, set, ave), set, con);
                                }
                            }

                            GADGET_CHECK_RETURN_FALSE(this->performPSIR(inputContainer, gmap, magIRContainer, magIRNoSCCContainer, PSIRContainer, psIRNoSCCContainer, magPDContainer, windowCenter, windowWidth));

                            if (windowCenter != -1)
                            {
                                wCenter[slc] = windowCenter;
                                wWidth[slc] = windowWidth;
                                psir_reconer_.compute_PSIR_windowing_ = false; // windowing are only computed once
                            }

                            for (con = 0; con < CON; con++)
                            {
                                magIR(cha, slc, con, phs, rep, 0, ave) = magIRContainer(0, con);
                                magIR(cha, slc, con, phs, rep, 0, ave).header_ = input(cha, slc, con, phs, rep, 0, ave).header_;
                                magIR(cha, slc, con, phs, rep, 0, ave).attrib_ = input(cha, slc, con, phs, rep, 0, ave).attrib_;
                                magIR(cha, slc, con, phs, rep, 0, ave).attrib_.set(GADGETRON_IMAGE_SCALE_RATIO, psir_reconer_.scale_factor_after_SCC_);
                                if (!TI_.empty()) magIR(cha, slc, con, phs, rep, 0, ave).attrib_.set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                                magIRNoSCC(cha, slc, con, phs, rep, 0, ave) = magIRNoSCCContainer(0, con);
                                magIRNoSCC(cha, slc, con, phs, rep, 0, ave).header_ = input(cha, slc, con, phs, rep, 0, ave).header_;
                                magIRNoSCC(cha, slc, con, phs, rep, 0, ave).attrib_ = input(cha, slc, con, phs, rep, 0, ave).attrib_;
                                if (!TI_.empty()) magIRNoSCC(cha, slc, con, phs, rep, 0, ave).attrib_.set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                                psIR(cha, slc, con, phs, rep, 0, ave) = PSIRContainer(0, con);
                                psIR(cha, slc, con, phs, rep, 0, ave).header_ = input(cha, slc, con, phs, rep, 0, ave).header_;
                                psIR(cha, slc, con, phs, rep, 0, ave).attrib_ = input(cha, slc, con, phs, rep, 0, ave).attrib_;
                                psIR(cha, slc, con, phs, rep, 0, ave).attrib_.set(GADGETRON_IMAGE_SCALE_RATIO, psir_reconer_.scale_factor_after_SCC_);
                                psIR(cha, slc, con, phs, rep, 0, ave).attrib_.set(GADGETRON_IMAGE_SCALE_OFFSET, psir_reconer_.offset_after_SCC_);
                                if (!TI_.empty()) psIR(cha, slc, con, phs, rep, 0, ave).attrib_.set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                                psIRNoSCC(cha, slc, con, phs, rep, 0, ave) = psIRNoSCCContainer(0, con);
                                psIRNoSCC(cha, slc, con, phs, rep, 0, ave).header_ = input(cha, slc, con, phs, rep, 0, ave).header_;
                                psIRNoSCC(cha, slc, con, phs, rep, 0, ave).attrib_ = input(cha, slc, con, phs, rep, 0, ave).attrib_;
                                if (!TI_.empty()) psIRNoSCC(cha, slc, con, phs, rep, 0, ave).attrib_.set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                                magPD(cha, slc, con, phs, rep, 0, ave) = magPDContainer(0, con);
                                magPD(cha, slc, con, phs, rep, 0, ave).header_ = input(cha, slc, con, phs, rep, 0, ave).header_;
                                magPD(cha, slc, con, phs, rep, 0, ave).attrib_ = input(cha, slc, con, phs, rep, 0, ave).attrib_;
                            }

                            magIRContainer.delete_data_on_destruct(false);
                            magIRNoSCCContainer.delete_data_on_destruct(false);
                            PSIRContainer.delete_data_on_destruct(false);
                            psIRNoSCCContainer.delete_data_on_destruct(false);
                            magPDContainer.delete_data_on_destruct(false);
                        }
                    }
                }
            }
        }

        return true;
    }

    bool GtImageMoCoAvePSIRGadget::performPSIR(ImageContainer2DType& input, ImageType& gmap, ImageContainer2DType& magIR, ImageContainer2DType& magIRNoSCC, ImageContainer2DType& psIR, ImageContainer2DType& psIRNoSCC, ImageContainer2DType& magPD, float& windowCenter, float& windowWidth)
    {
        try
        {
            magIR.clear();
            magIRNoSCC.clear();
            psIR.clear();
            psIRNoSCC.clear();
            magPD.clear();

            if (gmap.dimensions_equal(input(0, 0)))
            {
                psir_reconer_.gmap_ = gmap;
            }
            else
            {
                psir_reconer_.gmap_.clear();
            }

            GADGET_CHECK_RETURN_FALSE(psir_reconer_.performPSIRRecon(input));

            GADGET_CHECK_RETURN_FALSE(magIR.copyFrom(psir_reconer_.magIR_));
            GADGET_CHECK_RETURN_FALSE(magIRNoSCC.copyFrom(psir_reconer_.magIR_without_scc_));
            GADGET_CHECK_RETURN_FALSE(psIR.copyFrom(psir_reconer_.PSIR_));
            GADGET_CHECK_RETURN_FALSE(psIRNoSCC.copyFrom(psir_reconer_.PSIR_without_scc_));
            GADGET_CHECK_RETURN_FALSE(magPD.copyFrom(psir_reconer_.magPD_));

            windowCenter = psir_reconer_.window_center_;
            windowWidth = psir_reconer_.window_width_;
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GtImageMoCoAvePSIRGadget::performPSIR(ImageContainer2DType& input, ImageContainer2DType& magIR, ImageContainer2DType& psIR, ImageContainer2DType& magPD, float& windowCenter, float& windowWidth) ... ");
            return false;
        }

        return true;
    }

    int GtImageMoCoAvePSIRGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GtImageMoCoAvePSIRGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GtImageMoCoAvePSIRGadget)

}
