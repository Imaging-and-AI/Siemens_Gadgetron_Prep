
#include "GtImageMoCoAveGadget.h"
#include "mri_core_utility.h"

namespace Gadgetron { 

    GtImageMoCoAveGadget::GtImageMoCoAveGadget() : BaseClass()
    {
        moco_dim_ = DIM_Repetition;
        moco_cross_row_dim_ = DIM_Set;

        moco_ave_ = true;

        moco_cross_row_ = false;
        cross_row_same_reference_ = false;

        ref_moco_cross_row_ = 0;
        row_ref_pick_strategy_ = "Deformation";
        moco_quality_determined_in_ref_picking_ = true;

        strategy_ = GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE;
        dissimilarity_ = GT_IMAGE_DISSIMILARITY_LocalCCR;

        level_ = 3;
        iter_.resize(level_, 32);
        iter_[0] = 16;

        regularization_hilbert_strength_ = 12.0f;
        bidirectional_moco_ = true;

        level_cross_row_ = 3;
        iter_cross_row_.resize(level_cross_row_, 32);
        iter_cross_row_[0] = 16;

        dissimilarity_cross_row_ = GT_IMAGE_DISSIMILARITY_LocalCCR;
        regularization_hilbert_strength_cross_row_ = 24.0f;
        bidirectional_moco_cross_row_ = true;
        divergence_free_constraint_cross_row_ = false;

        dissimilarity_thres_ = 1e-5f;
        div_num_ = 2;
        inverse_deform_enforce_iter_ = 10;
        inverse_deform_enforce_weight_ = 0.5f;

        percentage_kept_for_averaging_ = 0.5;
        soft_averaging_ = true;

        send_ori_ = false;
        send_moco_ = false;
        send_moco_ave_ = true;

        moco_ave_keep_origial_image_number_ = false;
    }

    GtImageMoCoAveGadget::~GtImageMoCoAveGadget()
    {
    }

    bool GtImageMoCoAveGadget::readParameters()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(BaseClass::readParameters());

            GDEBUG_CONDITION_STREAM(this->verbose.value(), "------> GtImageMoCoAveGadget parameters <------");

            debugFolder_register_ = debugFolder_register.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "debugFolder_register_ is " << debugFolder_register_);

            if (!debugFolder_register_.empty())
            {
                Gadgetron::get_debug_folder_path(debugFolder_register_, debugFolder_register_fullPath_);
            }
            else
            {
                GDEBUG_STREAM("GtPlusImageMoCoAvePSIRGadget, debugFolder_register is not set ...");
            }

            verboseModeMOCO_ = verboseModeMOCO.value();
            GDEBUG_CONDITION_STREAM(verboseModeMOCO_, "verboseModeMOCO_ is " << verboseModeMOCO_);

            moco_dim_ = Gadgetron::get_ismrmrd_dim_from_name(moco_dim.value());
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "moco_dim_ is " << moco_dim.value());

            moco_cross_row_dim_ = Gadgetron::get_ismrmrd_dim_from_name(moco_cross_row_dim.value());
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "moco_cross_row_dim_ is " << moco_cross_row_dim.value());

            moco_ave_ = moco_ave.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "moco_ave_ is " << moco_ave_);

            moco_cross_row_ = moco_cross_row.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "moco_cross_row_ is " << moco_cross_row_);

            cross_row_same_reference_ = cross_row_same_reference.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "cross_row_same_reference_ is " << cross_row_same_reference_);

            ref_moco_cross_row_ = ref_moco_cross_row.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "ref_moco_cross_row_ is " << ref_moco_cross_row_);

            row_ref_pick_strategy_ = row_ref_pick_strategy.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "row_ref_pick_strategy_ is " << row_ref_pick_strategy_);

            moco_quality_determined_in_ref_picking_ = moco_quality_determined_in_ref_picking.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "moco_quality_determined_in_ref_picking_ is " << moco_quality_determined_in_ref_picking_);

            send_ori_ = send_ori.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "send_ori_ is " << send_ori_);

            send_moco_ = send_moco.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "send_moco_ is " << send_moco_);

            send_moco_ave_ = send_moco_ave.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "send_moco_ave_ is " << send_moco_ave_);

            moco_ave_keep_origial_image_number_ = moco_ave_keep_origial_image_number.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "moco_ave_keep_origial_image_number_ is " << moco_ave_keep_origial_image_number_);

            // -------------------------------------------------------

            strategy_ = Gadgetron::getImageRegContainerModeType(strategy.value());
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "strategy_ is " << strategy.value());

            dissimilarity_ = Gadgetron::getDissimilarityType(dissimilarity.value());
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "dissimilarity_ is " << dissimilarity.value());

            level_ = level.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "level_ is " << level_);

            iter_.resize(level_, 32);

            unsigned int ii;
            for (ii = 0; ii < level_; ii++)
            {
                if (ii == 0)
                    iter_[ii] = iter_0.value();
                else if (ii == 1)
                    iter_[ii] = iter_1.value();
                else if (ii == 2)
                    iter_[ii] = iter_2.value();
                else if (ii == 3)
                    iter_[ii] = iter_3.value();
                else if (ii == 4)
                    iter_[ii] = iter_4.value();
                else if (ii == 5)
                    iter_[ii] = iter_5.value();
                else if (ii == 6)
                    iter_[ii] = iter_6.value();
                else if (ii == 7)
                    iter_[ii] = iter_7.value();
                else if (ii == 8)
                    iter_[ii] = iter_8.value();
                else if (ii == 9)
                    iter_[ii] = iter_9.value();
                else if (ii == 10)
                    iter_[ii] = iter_10.value();
                else if (ii == 11)
                    iter_[ii] = iter_11.value();
                else if (ii == 12)
                    iter_[ii] = iter_12.value();
                else if (ii == 13)
                    iter_[ii] = iter_13.value();
                else if (ii == 14)
                    iter_[ii] = iter_14.value();
                else if (ii == 15)
                    iter_[ii] = iter_15.value();

                GDEBUG_CONDITION_STREAM(this->verbose.value(), "iter_ for level " << ii << " is " << iter_[ii]);
            }

            regularization_hilbert_strength_ = regularization_hilbert_strength.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "regularization_hilbert_strength_ is " << regularization_hilbert_strength_);

            bidirectional_moco_ = bidirectional_moco.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "bidirectional_moco_ is " << bidirectional_moco_);

            // -------------------------------------------------------

            dissimilarity_cross_row_ = Gadgetron::getDissimilarityType(dissimilarity_cross_row.value());
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "dissimilarity_cross_row_ is " << dissimilarity_cross_row.value());

            level_cross_row_ = level_cross_row.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "level_cross_row_ is " << level_cross_row_);

            iter_cross_row_.resize(level_cross_row_, 32);

            for (ii = 0; ii < level_cross_row_; ii++)
            {
                if (ii == 0)
                    iter_cross_row_[ii] = iter_cross_row_0.value();
                else if (ii == 1)
                    iter_cross_row_[ii] = iter_cross_row_1.value();
                else if (ii == 2)
                    iter_cross_row_[ii] = iter_cross_row_2.value();
                else if (ii == 3)
                    iter_cross_row_[ii] = iter_cross_row_3.value();
                else if (ii == 4)
                    iter_cross_row_[ii] = iter_cross_row_4.value();
                else if (ii == 5)
                    iter_cross_row_[ii] = iter_cross_row_5.value();
                else if (ii == 6)
                    iter_cross_row_[ii] = iter_cross_row_6.value();
                else if (ii == 7)
                    iter_cross_row_[ii] = iter_cross_row_7.value();
                else if (ii == 8)
                    iter_cross_row_[ii] = iter_cross_row_8.value();
                else if (ii == 9)
                    iter_cross_row_[ii] = iter_cross_row_9.value();
                else if (ii == 10)
                    iter_cross_row_[ii] = iter_cross_row_10.value();
                else if (ii == 11)
                    iter_cross_row_[ii] = iter_cross_row_11.value();
                else if (ii == 12)
                    iter_cross_row_[ii] = iter_cross_row_12.value();
                else if (ii == 13)
                    iter_cross_row_[ii] = iter_cross_row_13.value();
                else if (ii == 14)
                    iter_cross_row_[ii] = iter_cross_row_14.value();
                else if (ii == 15)
                    iter_cross_row_[ii] = iter_cross_row_15.value();

                GDEBUG_CONDITION_STREAM(this->verbose.value(), "iter_cross_row_ for level " << ii << " is " << iter_cross_row_[ii]);
            }

            regularization_hilbert_strength_cross_row_ = regularization_hilbert_strength_cross_row.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "regularization_hilbert_strength_cross_row_ is " << regularization_hilbert_strength_cross_row_);

            bidirectional_moco_cross_row_ = bidirectional_moco_cross_row.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "bidirectional_moco_cross_row_ is " << bidirectional_moco_cross_row_);

            divergence_free_constraint_cross_row_ = divergence_free_constraint_cross_row.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "divergence_free_constraint_cross_row_ is " << divergence_free_constraint_cross_row_);

            // -------------------------------------------------------

            dissimilarity_thres_ = dissimilarity_thres.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "dissimilarity_thres_ is " << dissimilarity_thres_);

            div_num_ = div_num.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "div_num_ is " << div_num_);

            inverse_deform_enforce_iter_ = inverse_deform_enforce_iter.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "inverse_deform_enforce_iter_ is " << inverse_deform_enforce_iter_);

            inverse_deform_enforce_weight_ = inverse_deform_enforce_weight.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "inverse_deform_enforce_weight_ is " << inverse_deform_enforce_weight_);

            GDEBUG_CONDITION_STREAM(this->verbose.value(), "-----------------------------------------------");

            percentage_kept_for_averaging_ = percentage_kept_for_averaging.value();
            if (percentage_kept_for_averaging_ == 0) percentage_kept_for_averaging_ = 0.5;
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "percentage_kept_for_averaging_ is " << percentage_kept_for_averaging_);

            soft_averaging_ = soft_averaging.value();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "soft_averaging_ is " << soft_averaging_);

            GDEBUG_CONDITION_STREAM(this->verbose.value(), "-----------------------------------------------");
        }
        catch (...)
        {
            GERROR_STREAM("Errors in GtImageMoCoAveGadget::readParameters() ... ");
            return false;
        }

        return true;
    }

    int GtImageMoCoAveGadget::process_config(ACE_Message_Block* mb)
    {
        // read in parameters from the xml
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);
        GADGET_CHECK_RETURN(this->readParameters(), GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...) {
            GDEBUG("Error parsing ISMRMRD Header");
            throw;
            return GADGET_FAIL;
        }

        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
        size_t max_segment = e_limits.segment ? e_limits.segment->maximum : 0;

        // ---------------------------------------------------------
        // if it is a semgneted acquisition, turn off the cross-row moco
        if (max_segment>0)
        {
            GDEBUG_STREAM("Turn off Cross-Row moco, a segmented acquisition is detected ... ");
            moco_cross_row_ = false;

        }
        // ---------------------------------------------------------

        size_t ii;

        // ---------------------------------------------------------
        // set the moco 2D
        // ---------------------------------------------------------
        mocoer2D_.debugFolder_ = this->debug_folder_full_path_;
        mocoer2D_.verbose_ = this->verbose.value();
        mocoer2D_.performTiming_ = this->perform_timing.value();

        mocoer2D_.cross_row_reg_ = moco_cross_row_;
        mocoer2D_.cross_row_reference_ = cross_row_same_reference_;

        mocoer2D_.row_ref_pick_strategy_ = row_ref_pick_strategy_;
        mocoer2D_.moco_quality_determined_in_ref_picking_ = moco_quality_determined_in_ref_picking_;

        mocoer2D_.percentage_kept_for_averaging_ = percentage_kept_for_averaging_;
        mocoer2D_.soft_averaging_ = soft_averaging_;

        mocoer2D_.register_.setDefaultParameters(level_, false);

        mocoer2D_.register_.container_reg_mode_ = GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE;
        mocoer2D_.register_.bg_value_ = -1;
        mocoer2D_.register_.debugFolder_ = this->debugFolder_register_fullPath_;
        mocoer2D_.register_.verbose_ = this->verboseModeMOCO_;
        mocoer2D_.register_.container_reg_transformation_ = (bidirectional_moco_ ? GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL : GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD);
        mocoer2D_.register_.max_iter_num_pyramid_level_ = iter_;

        mocoer2D_.register_.regularization_hilbert_strength_pyramid_level_.clear();
        mocoer2D_.register_.regularization_hilbert_strength_pyramid_level_.resize(level_);
        for (ii = 0; ii<level_; ii++)
        {
            mocoer2D_.register_.regularization_hilbert_strength_pyramid_level_[ii].resize(2, regularization_hilbert_strength_);
        }

        mocoer2D_.register_.dissimilarity_type_ = dissimilarity_;
        mocoer2D_.register_.dissimilarity_thres_pyramid_level_.clear();
        mocoer2D_.register_.dissimilarity_thres_pyramid_level_.resize(level_, dissimilarity_thres_);

        mocoer2D_.register_.inverse_deform_enforce_iter_pyramid_level_.clear();
        mocoer2D_.register_.inverse_deform_enforce_iter_pyramid_level_.resize(level_, inverse_deform_enforce_iter_);

        mocoer2D_.register_.inverse_deform_enforce_weight_pyramid_level_.clear();
        mocoer2D_.register_.inverse_deform_enforce_weight_pyramid_level_.resize(level_, inverse_deform_enforce_weight_);

        mocoer2D_.register_.div_num_pyramid_level_.clear();
        mocoer2D_.register_.div_num_pyramid_level_.resize(level_, div_num_);

        // ---------------------------------------------------------

        mocoer2D_.register_cross_row_.setDefaultParameters(level_cross_row_, false);
        mocoer2D_.register_cross_row_.container_reg_mode_ = GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE;
        mocoer2D_.register_cross_row_.bg_value_ = -1;
        mocoer2D_.register_cross_row_.debugFolder_ = this->debugFolder_register_fullPath_;
        mocoer2D_.register_cross_row_.verbose_ = this->verboseModeMOCO_;
        mocoer2D_.register_cross_row_.container_reg_transformation_ = (bidirectional_moco_cross_row_ ? GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL : GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD);
        mocoer2D_.register_cross_row_.max_iter_num_pyramid_level_ = iter_cross_row_;
        mocoer2D_.register_cross_row_.apply_divergence_free_constraint_ = divergence_free_constraint_cross_row_;

        mocoer2D_.register_cross_row_.regularization_hilbert_strength_pyramid_level_.clear();
        mocoer2D_.register_cross_row_.regularization_hilbert_strength_pyramid_level_.resize(level_cross_row_);
        for (ii = 0; ii<level_cross_row_; ii++)
        {
            mocoer2D_.register_cross_row_.regularization_hilbert_strength_pyramid_level_[ii].resize(2, regularization_hilbert_strength_cross_row_);
        }

        mocoer2D_.register_cross_row_.dissimilarity_type_ = dissimilarity_cross_row_;
        mocoer2D_.register_cross_row_.dissimilarity_thres_pyramid_level_.clear();
        mocoer2D_.register_cross_row_.dissimilarity_thres_pyramid_level_.resize(level_cross_row_, dissimilarity_thres_);

        mocoer2D_.register_cross_row_.inverse_deform_enforce_iter_pyramid_level_.clear();
        mocoer2D_.register_cross_row_.inverse_deform_enforce_iter_pyramid_level_.resize(level_cross_row_, inverse_deform_enforce_iter_);

        mocoer2D_.register_cross_row_.inverse_deform_enforce_weight_pyramid_level_.clear();
        mocoer2D_.register_cross_row_.inverse_deform_enforce_weight_pyramid_level_.resize(level_cross_row_, inverse_deform_enforce_weight_);

        mocoer2D_.register_cross_row_.div_num_pyramid_level_.clear();
        mocoer2D_.register_cross_row_.div_num_pyramid_level_.resize(level_cross_row_, div_num_);

        // set up the ref selection mocoer from cross-row mocoer
        mocoer2D_.register_cross_row_ref_selection_ = mocoer2D_.register_cross_row_;
        mocoer2D_.register_cross_row_ref_selection_.max_iter_num_pyramid_level_[0] = 1;
        mocoer2D_.register_cross_row_ref_selection_.container_reg_transformation_ = GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD;

        return GADGET_OK;
    }

    int GtImageMoCoAveGadget::processImageBuffer(ImageBufferType& ori)
    {
        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::processImageBuffer(...) starts ... ");

        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        std::vector<size_t> dims = ori.get_dimensions();
        GDEBUG_CONDITION_STREAM(this->verbose.value(), "[Cha Slice Con Phase Rep Set Ave] = [" << dims[0] << " " << dims[1] << " " << dims[2] << " "
            << dims[3] << " " << dims[4] << " " << dims[5] << " "
            << dims[6] << "]");

        if (send_ori_)
        {
            GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_, processStr, dataRole), GADGET_FAIL);
        }

        ImageBufferType moco, ave;

        if (moco_ave_)
        {
            if (this->perform_timing.value()) { gt_timer_.start("Perform MOCO AVE ... "); }
            bool mocoPerformed(true);
            GADGET_CHECK_RETURN(this->performMoCoAveraging(ori, moco, ave, mocoPerformed), GADGET_FAIL);
            if (this->perform_timing.value()) { gt_timer_.stop(); }

            // send out the images
            if (send_moco_)
            {
                if (mocoPerformed)
                {
                    processStr.push_back(GADGETRON_IMAGE_MOCO);
                    dataRole.push_back(GADGETRON_IMAGE_MOCO);
                }

                GADGET_CHECK_RETURN(this->sendOutImages(moco, image_series_num_ + 1, processStr, dataRole), GADGET_FAIL);
            }

            if (send_moco_ave_)
            {
                processStr.clear();
                if (mocoPerformed) processStr.push_back(GADGETRON_IMAGE_MOCO);
                processStr.push_back(GADGETRON_IMAGE_AVE);

                dataRole.clear();
                if (mocoPerformed) dataRole.push_back(GADGETRON_IMAGE_MOCO);
                dataRole.push_back(GADGETRON_IMAGE_AVE);
                GADGET_CHECK_RETURN(this->sendOutImages(ave, image_series_num_ + 2, processStr, dataRole), GADGET_FAIL);
            }

            // release the moco images
            GADGET_CHECK_RETURN(this->releaseImageBuffer(moco), GADGET_FAIL);

            // release the moco ave images
            GADGET_CHECK_RETURN(this->releaseImageBuffer(ave), GADGET_FAIL);
        }

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::process(...) ends ... ");

        return GADGET_OK;
    }

    bool GtImageMoCoAveGadget::performMoCoAveraging(ImageBufferType& ori, ImageBufferType& moco, ImageBufferType& aveRes, bool& mocoPerformed)
    {
        size_t CHA = ori.get_size(0);
        size_t SLC = ori.get_size(1);
        size_t CON = ori.get_size(2);
        size_t PHS = ori.get_size(3);
        size_t REP = ori.get_size(4);
        size_t SET = ori.get_size(5);
        size_t AVE = ori.get_size(6);

        std::vector<size_t> dim;
        ori.get_dimensions(dim);

        moco.create(dim);
        mocoPerformed = true;

        if (moco_dim_ == DIM_Repetition)
        {
            dim[4] = 1;
            if (REP == 1) mocoPerformed = false;
        }
        else if (moco_dim_ == DIM_Phase)
        {
            dim[3] = 1;
            if (PHS == 1) mocoPerformed = false;
        }
        else if (moco_dim_ == DIM_Average)
        {
            dim[6] = 1;
            if (AVE == 1) mocoPerformed = false;
        }
        aveRes.create(dim);

        size_t cha, slc, con, phs, rep, set, ave;

        for (cha = 0; cha<CHA; cha++)
        {
            if (moco_dim_ == DIM_Repetition)
            {
                if (moco_cross_row_dim_ == DIM_Set)
                {
                    ImageContainer2DType inputContainer, mocoContainer, aveContainer;
                    std::vector<size_t> cols(SET*SLC, REP);
                    inputContainer.create(cols, false);

                    for (ave = 0; ave<AVE; ave++)
                    {
                        for (phs = 0; phs<PHS; phs++)
                        {
                            // assemble image container for the first contrast
                            for (con = 0; con<CON; con++)
                            {
                                for (slc = 0; slc<SLC; slc++)
                                {
                                    // fill the image container
                                    for (set = 0; set<SET; set++)
                                    {
                                        for (rep = 0; rep<REP; rep++)
                                        {
                                            inputContainer.set(&ori(cha, slc, con, phs, rep, set, ave), set + slc * SET, rep);
                                        }
                                    }
                                }

                                if (con == 0)
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for first contrast ... ");
                                    this->performMoCoAveraging2D(inputContainer, mocoContainer, aveContainer, SET, SLC);
                                }
                                else
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for other contrasts ... ");
                                    this->applyMoCoAveraging2D(inputContainer, mocoContainer, aveContainer);
                                }

                                mocoContainer.delete_data_on_destruct(false);
                                aveContainer.delete_data_on_destruct(false);

                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (set = 0; set<SET; set++)
                                    {
                                        for (rep = 0; rep<REP; rep++)
                                        {
                                            moco(cha, slc, con, phs, rep, set, ave) = mocoContainer(set + slc * SET, rep);
                                            moco(cha, slc, con, phs, rep, set, ave).header_ = ori(cha, slc, con, phs, rep, set, ave).header_;
                                            moco(cha, slc, con, phs, rep, set, ave).attrib_ = ori(cha, slc, con, phs, rep, set, ave).attrib_;
                                        }

                                        aveRes(cha, slc, con, phs, 0, set, ave) = aveContainer(0, set + slc * SET);
                                        aveRes(cha, slc, con, phs, 0, set, ave).header_ = ori(cha, slc, con, phs, 0, set, ave).header_;
                                        aveRes(cha, slc, con, phs, 0, set, ave).attrib_ = ori(cha, slc, con, phs, 0, set, ave).attrib_;

                                        if (!moco_ave_keep_origial_image_number_)
                                        {
                                            // aveRes(cha, slc, con, phs, 0, set, ave).attrib_.set(ISMRMRD_IMAGE_repetition, 0L);
                                            aveRes(cha, slc, con, phs, 0, set, ave).header_.repetition = 0;
                                            aveRes(cha, slc, con, phs, 0, set, ave).attrib_.set(GADGETRON_IMAGENUMBER, 0L);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else if (moco_cross_row_dim_ == DIM_Contrast)
                {
                    ImageContainer2DType inputContainer, mocoContainer, aveContainer;
                    std::vector<size_t> cols(CON*SLC, REP);
                    inputContainer.create(cols, false);

                    for (ave = 0; ave<AVE; ave++)
                    {
                        for (phs = 0; phs<PHS; phs++)
                        {
                            // assemble image container for the first contrast
                            for (set = 0; set<SET; set++)
                            {
                                // fill the image container
                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (con = 0; con<CON; con++)
                                    {
                                        for (rep = 0; rep<REP; rep++)
                                        {
                                            inputContainer.set(&ori(cha, slc, con, phs, rep, set, ave), con + slc * CON, rep);
                                        }
                                    }
                                }

                                if (set == 0)
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for first set ... ");
                                    this->performMoCoAveraging2D(inputContainer, mocoContainer, aveContainer, CON, SLC);
                                }
                                else
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::apply motion correction for other sets ... ");
                                    this->applyMoCoAveraging2D(inputContainer, mocoContainer, aveContainer);
                                }

                                mocoContainer.delete_data_on_destruct(false);
                                aveContainer.delete_data_on_destruct(false);

                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (con = 0; con<CON; con++)
                                    {
                                        for (rep = 0; rep<REP; rep++)
                                        {
                                            moco(cha, slc, con, phs, rep, set, ave) = mocoContainer(con, rep);
                                            moco(cha, slc, con, phs, rep, set, ave).header_ = ori(cha, slc, con, phs, rep, set, ave).header_;
                                            moco(cha, slc, con, phs, rep, set, ave).attrib_ = ori(cha, slc, con, phs, rep, set, ave).attrib_;
                                        }

                                        aveRes(cha, slc, con, phs, 0, set, ave) = aveContainer(0, con + slc * CON);
                                        aveRes(cha, slc, con, phs, 0, set, ave).header_ = ori(cha, slc, con, phs, 0, set, ave).header_;
                                        aveRes(cha, slc, con, phs, 0, set, ave).attrib_ = ori(cha, slc, con, phs, 0, set, ave).attrib_;

                                        if (!moco_ave_keep_origial_image_number_)
                                        {
                                            // aveRes(cha, slc, con, phs, 0, set, ave).attrib_.set(ISMRMRD_IMAGE_repetition, 0L);
                                            aveRes(cha, slc, con, phs, 0, set, ave).header_.repetition = 0;
                                            aveRes(cha, slc, con, phs, 0, set, ave).attrib_.set(GADGETRON_IMAGENUMBER, 0L);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (moco_dim_ == DIM_Phase)
            {
                if (moco_cross_row_dim_ == DIM_Set)
                {
                    ImageContainer2DType inputContainer, mocoContainer, aveContainer;
                    std::vector<size_t> cols(SET*SLC, PHS);
                    inputContainer.create(cols, false);

                    for (ave = 0; ave<AVE; ave++)
                    {
                        for (rep = 0; rep<REP; rep++)
                        {
                            // assemble image container for the first contrast
                            for (con = 0; con<CON; con++)
                            {
                                // fill the image container
                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (set = 0; set<SET; set++)
                                    {
                                        for (phs = 0; phs<PHS; phs++)
                                        {
                                            inputContainer.set(&ori(cha, slc, con, phs, rep, set, ave), set + slc * SET, phs);
                                        }
                                    }
                                }

                                if (con == 0)
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for first contrast ... ");
                                    this->performMoCoAveraging2D(inputContainer, mocoContainer, aveContainer, SET, SLC);
                                }
                                else
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::apply motion correction for other contrasts ... ");
                                    this->applyMoCoAveraging2D(inputContainer, mocoContainer, aveContainer);
                                }

                                mocoContainer.delete_data_on_destruct(false);
                                aveContainer.delete_data_on_destruct(false);

                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (set = 0; set<SET; set++)
                                    {
                                        for (phs = 0; phs<PHS; phs++)
                                        {
                                            moco(cha, slc, con, phs, rep, set, ave) = mocoContainer(set, phs);
                                            moco(cha, slc, con, phs, rep, set, ave).header_ = ori(cha, slc, con, phs, rep, set, ave).header_;
                                            moco(cha, slc, con, phs, rep, set, ave).attrib_ = ori(cha, slc, con, phs, rep, set, ave).attrib_;
                                        }

                                        aveRes(cha, slc, con, 0, rep, set, ave) = aveContainer(0, set + slc * SET);
                                        aveRes(cha, slc, con, 0, rep, set, ave).header_ = ori(cha, slc, con, 0, rep, set, ave).header_;
                                        aveRes(cha, slc, con, 0, rep, set, ave).attrib_ = ori(cha, slc, con, 0, rep, set, ave).attrib_;

                                        if (!moco_ave_keep_origial_image_number_)
                                        {
                                            // aveRes(cha, slc, con, 0, rep, set, ave).attrib_.set(ISMRMRD_IMAGE_phase, 0L);
                                            aveRes(cha, slc, con, 0, rep, set, ave).header_.phase = 0;
                                            aveRes(cha, slc, con, 0, rep, set, ave).attrib_.set(GADGETRON_IMAGENUMBER, 0L);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else if (moco_cross_row_dim_ == DIM_Contrast)
                {
                    ImageContainer2DType inputContainer, mocoContainer, aveContainer;
                    std::vector<size_t> cols(CON*SLC, PHS);
                    inputContainer.create(cols, false);

                    for (ave = 0; ave<AVE; ave++)
                    {
                        for (rep = 0; rep<REP; rep++)
                        {
                            // assemble image container for the first contrast
                            for (set = 0; set<SET; set++)
                            {
                                // fill the image container
                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (con = 0; con<CON; con++)
                                    {
                                        for (phs = 0; phs<PHS; phs++)
                                        {
                                            inputContainer.set(&ori(cha, slc, con, phs, rep, set, ave), con + slc * CON, phs);
                                        }
                                    }
                                }

                                if (set == 0)
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for first set ... ");
                                    this->performMoCoAveraging2D(inputContainer, mocoContainer, aveContainer, CON, SLC);
                                }
                                else
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::apply motion correction for other sets ... ");
                                    this->applyMoCoAveraging2D(inputContainer, mocoContainer, aveContainer);
                                }

                                mocoContainer.delete_data_on_destruct(false);
                                aveContainer.delete_data_on_destruct(false);

                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (con = 0; con<CON; con++)
                                    {
                                        for (phs = 0; phs<PHS; phs++)
                                        {
                                            moco(cha, slc, con, phs, rep, set, ave) = mocoContainer(con, phs);
                                            moco(cha, slc, con, phs, rep, set, ave).header_ = ori(cha, slc, con, phs, rep, set, ave).header_;
                                            moco(cha, slc, con, phs, rep, set, ave).attrib_ = ori(cha, slc, con, phs, rep, set, ave).attrib_;
                                        }

                                        aveRes(cha, slc, con, 0, rep, set, ave) = aveContainer(0, con + slc * CON);
                                        aveRes(cha, slc, con, 0, rep, set, ave).header_ = ori(cha, slc, con, 0, rep, set, ave).header_;
                                        aveRes(cha, slc, con, 0, rep, set, ave).attrib_ = ori(cha, slc, con, 0, rep, set, ave).attrib_;

                                        if (!moco_ave_keep_origial_image_number_)
                                        {
                                            // aveRes(cha, slc, con, 0, rep, set, ave).attrib_.set(ISMRMRD_IMAGE_phase, 0L);
                                            aveRes(cha, slc, con, 0, rep, set, ave).header_.phase = 0;
                                            aveRes(cha, slc, con, 0, rep, set, ave).attrib_.set(GADGETRON_IMAGENUMBER, 0L);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (moco_dim_ == DIM_Average)
            {
                if (moco_cross_row_dim_ == DIM_Set)
                {
                    ImageContainer2DType inputContainer, mocoContainer, aveContainer;
                    std::vector<size_t> cols(SET*SLC, AVE);
                    inputContainer.create(cols, false);

                    for (rep = 0; rep<REP; rep++)
                    {
                        for (phs = 0; phs<PHS; phs++)
                        {
                            // assemble image container for the first contrast
                            for (con = 0; con<CON; con++)
                            {
                                for (slc = 0; slc<SLC; slc++)
                                {
                                    // fill the image container
                                    for (set = 0; set<SET; set++)
                                    {
                                        for (ave = 0; ave<AVE; ave++)
                                        {
                                            inputContainer.set(&ori(cha, slc, con, phs, rep, set, ave), set + slc * SET, ave);
                                        }
                                    }
                                }

                                if (con == 0)
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for first contrast ... ");
                                    this->performMoCoAveraging2D(inputContainer, mocoContainer, aveContainer, SET, SLC);
                                }
                                else
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::apply motion correction for other contrasts ... ");
                                    this->applyMoCoAveraging2D(inputContainer, mocoContainer, aveContainer);
                                }

                                mocoContainer.delete_data_on_destruct(false);
                                aveContainer.delete_data_on_destruct(false);

                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (set = 0; set<SET; set++)
                                    {
                                        for (ave = 0; ave<AVE; ave++)
                                        {
                                            moco(cha, slc, con, phs, rep, set, ave) = mocoContainer(set + slc * SET, ave);
                                            moco(cha, slc, con, phs, rep, set, ave).header_ = ori(cha, slc, con, phs, rep, set, ave).header_;
                                            moco(cha, slc, con, phs, rep, set, ave).attrib_ = ori(cha, slc, con, phs, rep, set, ave).attrib_;
                                        }

                                        aveRes(cha, slc, con, phs, rep, set, 0) = aveContainer(0, set + slc * SET);
                                        aveRes(cha, slc, con, phs, rep, set, 0).header_ = ori(cha, slc, con, phs, rep, set, 0).header_;
                                        aveRes(cha, slc, con, phs, rep, set, 0).attrib_ = ori(cha, slc, con, phs, rep, set, 0).attrib_;

                                        if (!moco_ave_keep_origial_image_number_)
                                        {
                                            // aveRes(cha, slc, con, phs, rep, set, 0).attrib_.set(ISMRMRD_IMAGE_average, 0L);
                                            aveRes(cha, slc, con, phs, rep, set, 0).header_.average = 0;
                                            aveRes(cha, slc, con, phs, rep, set, 0).attrib_.set(GADGETRON_IMAGENUMBER, 0L);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else if (moco_cross_row_dim_ == DIM_Contrast)
                {
                    ImageContainer2DType inputContainer, mocoContainer, aveContainer;
                    std::vector<size_t> cols(CON*SLC, AVE);
                    inputContainer.create(cols, false);

                    for (rep = 0; rep<REP; rep++)
                    {
                        for (phs = 0; phs<PHS; phs++)
                        {
                            // assemble image container for the first contrast
                            for (set = 0; set<SET; set++)
                            {
                                // fill the image container
                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (con = 0; con<CON; con++)
                                    {
                                        for (ave = 0; ave<AVE; ave++)
                                        {
                                            inputContainer.set(&ori(cha, slc, con, phs, rep, set, ave), con + slc * CON, ave);
                                        }
                                    }
                                }

                                if (set == 0)
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for first set ... ");
                                    this->performMoCoAveraging2D(inputContainer, mocoContainer, aveContainer, CON, SLC);
                                }
                                else
                                {
                                    GDEBUG_CONDITION_STREAM(this->verbose.value(), "GtImageMoCoAveGadget::perform motion correction for other sets ... ");
                                    this->applyMoCoAveraging2D(inputContainer, mocoContainer, aveContainer);
                                }

                                mocoContainer.delete_data_on_destruct(false);
                                aveContainer.delete_data_on_destruct(false);

                                for (slc = 0; slc<SLC; slc++)
                                {
                                    for (con = 0; con<CON; con++)
                                    {
                                        for (ave = 0; ave<AVE; ave++)
                                        {
                                            moco(cha, slc, con, phs, rep, set, ave) = mocoContainer(con, ave);
                                            moco(cha, slc, con, phs, rep, set, ave).header_ = ori(cha, slc, con, phs, rep, set, ave).header_;
                                            moco(cha, slc, con, phs, rep, set, ave).attrib_ = ori(cha, slc, con, phs, rep, set, ave).attrib_;
                                        }

                                        aveRes(cha, slc, con, phs, rep, set, 0) = aveContainer(0, con + slc * CON);
                                        aveRes(cha, slc, con, phs, rep, set, 0).header_ = ori(cha, slc, con, phs, rep, set, 0).header_;
                                        aveRes(cha, slc, con, phs, rep, set, 0).attrib_ = ori(cha, slc, con, phs, rep, set, 0).attrib_;

                                        if (!moco_ave_keep_origial_image_number_)
                                        {
                                            // aveRes(cha, slc, con, phs, rep, set, 0).attrib_.set(ISMRMRD_IMAGE_average, 0L);
                                            aveRes(cha, slc, con, phs, rep, set, 0).header_.average = 0;
                                            aveRes(cha, slc, con, phs, rep, set, 0).attrib_.set(GADGETRON_IMAGENUMBER, 0L);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    bool GtImageMoCoAveGadget::performMoCoAveraging2D(ImageContainer2DType& inputContainer, ImageContainer2DType& mocoContainer, ImageContainer2DType& aveContainer, size_t groupSize, size_t groupNum)
    {
        try
        {
            size_t R = inputContainer.rows();
            std::vector<size_t> cols = inputContainer.cols();

            mocoContainer.clear();
            aveContainer.clear();

            if (!mag_container2D_.dimensions_equal_container(inputContainer))
            {
                mag_container2D_.clear();
                mag_container2D_.create(cols);
            }

            size_t r, c;
            for (r = 0; r<R; r++)
            {
                for (c = 0; c<cols[r]; c++)
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::abs(inputContainer(r, c), mag_container2D_(r, c)));
                }
            }

            this->exportImageContainer2D(mag_container2D_, "mag_for_moco");

            if (groupNum > 1)
            {
                GADGET_CHECK_RETURN_FALSE(groupSize*groupNum == R);

                mocoer2D_.row_group_index_.resize(R, 0);
                size_t ii, jj;
                for (jj = 0; jj<groupNum; jj++)
                {
                    for (ii = 0; ii<groupSize; ii++)
                    {
                        mocoer2D_.row_group_index_[ii + jj * groupSize] = jj;
                    }
                }
            }

            GADGET_CHECK_RETURN_FALSE(mocoer2D_.computeMoCoAveraging(mag_container2D_));

            this->exportImageContainer2D(mocoer2D_.register_.warped_container_, "mag_after_moco");
            this->exportImageContainer2D(mocoer2D_.averaged_, "mag_ave");

            if (moco_cross_row_)
            {
                if (groupNum <= 1)
                {
                    GADGET_CHECK_RETURN_FALSE(mocoer2D_.performCrossRowMoCo(ref_moco_cross_row_));
                }
                else
                {
                    std::vector<unsigned int> refRows(R);

                    size_t ii, jj;
                    for (jj = 0; jj<groupNum; jj++)
                    {
                        for (ii = 0; ii<groupSize; ii++)
                        {
                            refRows[ii + jj * groupSize] = (unsigned int)(ref_moco_cross_row_ + jj * groupSize);
                        }
                    }

                    GADGET_CHECK_RETURN_FALSE(mocoer2D_.performCrossRowMoCo(refRows));
                }

                this->exportImageContainer2D(mocoer2D_.averaged_, "mag_ave_cross_moco");
            }

            GADGET_CHECK_RETURN_FALSE(mocoContainer.copyFrom(inputContainer));
            GADGET_CHECK_RETURN_FALSE(mocoer2D_.applyMoCoAveraging(inputContainer, mocoContainer, aveContainer));

            this->exportImageContainer2D(mocoContainer, "moco");
            this->exportImageContainer2D(aveContainer, "ave");

            if (moco_cross_row_)
            {
                ImageContainer2DType averagedCxContainerMoCo;

                GADGET_CHECK_RETURN_FALSE(averagedCxContainerMoCo.copyFrom(aveContainer));

                GADGET_CHECK_RETURN_FALSE(mocoer2D_.applyCrossRowMoCo(aveContainer, averagedCxContainerMoCo));

                GADGET_CHECK_RETURN_FALSE(aveContainer.copyFrom(averagedCxContainerMoCo));

                this->exportImageContainer2D(aveContainer, "ave_cross_moco");
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GtImageMoCoAveGadget::performMoCoAveraging2D(...) ... ");
            return false;
        }

        return true;
    }

    bool GtImageMoCoAveGadget::applyMoCoAveraging2D(ImageContainer2DType& inputContainer, ImageContainer2DType& mocoContainer, ImageContainer2DType& aveContainer)
    {
        try
        {
            mocoContainer.clear();
            aveContainer.clear();

            GADGET_CHECK_RETURN_FALSE(mocoContainer.copyFrom(inputContainer));
            GADGET_CHECK_RETURN_FALSE(mocoer2D_.applyMoCoAveraging(inputContainer, mocoContainer, aveContainer));

            this->exportImageContainer2D(mocoContainer, "moco");
            this->exportImageContainer2D(aveContainer, "ave");

            if (moco_cross_row_)
            {
                ImageContainer2DType averagedCxContainerMoCo;

                GADGET_CHECK_RETURN_FALSE(averagedCxContainerMoCo.copyFrom(aveContainer));

                GADGET_CHECK_RETURN_FALSE(mocoer2D_.applyCrossRowMoCo(aveContainer, averagedCxContainerMoCo));

                GADGET_CHECK_RETURN_FALSE(aveContainer.copyFrom(averagedCxContainerMoCo));

                this->exportImageContainer2D(mocoContainer, "ave_cross_moco");
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GtImageMoCoAveGadget::applyMoCoAveraging2D(...) ... ");
            return false;
        }

        return true;
    }

    int GtImageMoCoAveGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GtImageMoCoAveGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GtImageMoCoAveGadget)

}
