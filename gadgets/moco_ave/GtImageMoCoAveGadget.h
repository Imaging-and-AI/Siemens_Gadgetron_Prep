/** 
\file   GtImageMoCoAveGadget.h
\brief  The image moco and averaging gadget
        The MoCo dimension can be REP or PHS, the cross-series dimension is SET
        The MoCo averaging results can be applied across CON
\author Hui Xue
*/

#pragma once

#include "gadgetron_siemens_moco_ave_export.h"
#include <GenericImageReconGadget.h>
#include "core_moco_averaging.h"

namespace Gadgetron {

    class EXPORTGTGADGET GtImageMoCoAveGadget : public GenericImageReconGadget
    {
    public:
        GADGET_DECLARE(GtImageMoCoAveGadget);

        typedef GenericImageReconGadget BaseClass;

        typedef BaseClass::ValueType ValueType;
        typedef BaseClass::Image2DType ImageType;
        typedef BaseClass::Image2DBufferType ImageBufferType;
        typedef BaseClass::ImgArrayType ImgArrayType;

        GtImageMoCoAveGadget();
        ~GtImageMoCoAveGadget();

        virtual int close(unsigned long flags);

        /// dimension to perform moco+ave
        /// can be REP or PHS
        Gadgetron::IsmrmrdDIM moco_dim_;

        /// dimension to perform cross-row moco
        /// can be SET or CON
        Gadgetron::IsmrmrdDIM moco_cross_row_dim_;

        /// whether to perform moco_ave
        bool moco_ave_;

        /// whether to perform cross-row moco
        bool moco_cross_row_;

        /// all rows have the same reference
        bool cross_row_same_reference_;

        /// which row as the reference
        unsigned int ref_moco_cross_row_;

        /// strategy to pick the reference for every row
        /// "SSD" or "Deformation"
        std::string row_ref_pick_strategy_;

        /// whether to determine moco quality in the ref picking phase
        /// if true, and row_ref_pick_strategy_ == "SSD", the SSD values will be used as moco quality
        /// if true, and row_ref_pick_strategy_ == "Deformation", the deformation field will be used as moco quality
        bool moco_quality_determined_in_ref_picking_;

        /// whether to send original images
        bool send_ori_;

        /// whether to send moco images
        bool send_moco_;

        /// whether to send averaged images
        bool send_moco_ave_;

        /// whether sent averaged images keep the original image number
        bool moco_ave_keep_origial_image_number_;

        // debug folder for register
        std::string debugFolder_register_;
        std::string debugFolder_register_fullPath_;

    protected:

        GADGET_PROPERTY(verboseModeMOCO, bool, "Whether to print out MOCO iterations", false);

        GADGET_PROPERTY(debugFolder_register, std::string, "If set, the debug output will be written out for register", "");

        GADGET_PROPERTY_LIMITS(moco_dim, std::string, "Dimension for motion correction", "DIM_Repetition",
            GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
            "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

        GADGET_PROPERTY_LIMITS(moco_cross_row_dim, std::string, "Dimension for cross-row motion correction", "DIM_Set",
            GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
            "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

        GADGET_PROPERTY(moco_ave, bool, "Whether to perform motion correction and averaging", true);
        GADGET_PROPERTY(moco_cross_row, bool, "Whether to perform motion correction cross-row", true);
        GADGET_PROPERTY(cross_row_same_reference, bool, "Whether all rows have the same reference frame", false);

        GADGET_PROPERTY(ref_moco_cross_row, int, "If cross-row MOCO is performed, which row is selected as the reference", 0);

        GADGET_PROPERTY_LIMITS(row_ref_pick_strategy, std::string, "Strategy to pick reference for rows", "Deformation",
            GadgetPropertyLimitsEnumeration, "SSD", "Deformation");

        GADGET_PROPERTY(moco_quality_determined_in_ref_picking, bool, "Whether to determine the moco quality while picking the reference", true);

        GADGET_PROPERTY(send_ori, bool, "Whether to send original images", false);
        GADGET_PROPERTY(send_moco, bool, "Whether to send moco images", false);
        GADGET_PROPERTY(send_moco_ave, bool, "Whether to send moco ave images", true);

        GADGET_PROPERTY(moco_ave_keep_origial_image_number, bool, "Whether to keep original image number for the moco ave images", false);

        // -------------------------------------------------------

        GADGET_PROPERTY_LIMITS(strategy, std::string, "MOCO strategy", "FixedReference",
            GadgetPropertyLimitsEnumeration, "FixedReference", "Progressive");

        GADGET_PROPERTY_LIMITS(dissimilarity, std::string, "Image dissimilarity measures", "LocalCCR",
            GadgetPropertyLimitsEnumeration, "SSD", "LocalCCR", "MutualInformation");

        GADGET_PROPERTY(level, int, "Number of levels for the multi-resolution pyramid for motion correction", 4);
        GADGET_PROPERTY(iter_0, int, "Number of moco iterations", 1);
        GADGET_PROPERTY(iter_1, int, "Number of moco iterations", 32);
        GADGET_PROPERTY(iter_2, int, "Number of moco iterations", 64);
        GADGET_PROPERTY(iter_3, int, "Number of moco iterations", 64);
        GADGET_PROPERTY(iter_4, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_5, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_6, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_7, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_8, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_9, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_10, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_11, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_12, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_13, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_14, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_15, int, "Number of moco iterations", 100);

        GADGET_PROPERTY(regularization_hilbert_strength, double, "Hilbert stregth for moco", 12.0);
        GADGET_PROPERTY(bidirectional_moco, bool, "Whether to apply bidirectional moco", false);

        // -------------------------------------------------------

        GADGET_PROPERTY_LIMITS(dissimilarity_cross_row, std::string, "Image dissimilarity measures for cross-row moco", "LocalCCR",
            GadgetPropertyLimitsEnumeration, "SSD", "LocalCCR", "MutualInformation");

        GADGET_PROPERTY(level_cross_row, int, "Number of levels for the multi-resolution pyramid for cross-row motion correction", 4);
        GADGET_PROPERTY(iter_cross_row_0, int, "Number of moco iterations", 16);
        GADGET_PROPERTY(iter_cross_row_1, int, "Number of moco iterations", 32);
        GADGET_PROPERTY(iter_cross_row_2, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_3, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_4, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_5, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_6, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_7, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_8, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_9, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_10, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_11, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_12, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_13, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_14, int, "Number of moco iterations", 100);
        GADGET_PROPERTY(iter_cross_row_15, int, "Number of moco iterations", 100);

        GADGET_PROPERTY(regularization_hilbert_strength_cross_row, double, "Hilbert stregth for cross-row moco", 36.0);
        GADGET_PROPERTY(bidirectional_moco_cross_row, bool, "Whether to apply bidirectional moco for cross-row", true);
        GADGET_PROPERTY(divergence_free_constraint_cross_row, bool, "Whether to apply divergence free constraint for cross-row moco", true);

        // -------------------------------------------------------

        GADGET_PROPERTY(dissimilarity_thres, double, "Threshold for image dissimilarity minimization", 0);
        GADGET_PROPERTY(div_num, int, "Number of sub-division search in minimization", 3);
        GADGET_PROPERTY(inverse_deform_enforce_iter, int, "For the bidirectional MOCO, the number of bidirectional iteration", 10);
        GADGET_PROPERTY(inverse_deform_enforce_weight, double, "For the bidirectional MOCO, the weight between forward and inverse MOCO", 0.5);

        GADGET_PROPERTY(percentage_kept_for_averaging, double, "Fraction of images kept for averaging", 0.5);
        GADGET_PROPERTY(soft_averaging, bool, "Whether to perform soft averaging", true);

        // -------------------------------------------------------

        // read in parameters
        bool readParameters();

        virtual int process_config(ACE_Message_Block* mb);
        virtual int processImageBuffer(ImageBufferType& ori);

        /// perform the moco+ave
        bool performMoCoAveraging(ImageBufferType& ori, ImageBufferType& moco, ImageBufferType& ave, bool& mocoPerformed);

        /// perform and apply moco+ave for 2D
        /// inputContainer has groupSize*groupNum rows
        /// if groupNum == 1, all rows in inputContainer are within one group
        bool performMoCoAveraging2D(ImageContainer2DType& inputContainer, ImageContainer2DType& mocoContainer, ImageContainer2DType& aveContainer, size_t groupSize, size_t groupNum);
        bool applyMoCoAveraging2D(ImageContainer2DType& inputContainer, ImageContainer2DType& mocoContainer, ImageContainer2DType& aveContainer);

        /// parameters for moco and ave
        /// moco mode
        Gadgetron::GT_IMAGE_REG_CONTAINER_MODE strategy_;
        /// image dissimilarity
        Gadgetron::GT_IMAGE_DISSIMILARITY dissimilarity_;
        /// number of resolution levels
        unsigned int level_;
        /// iterations for every level
        std::vector<unsigned int> iter_;
        /// regularization strength, hilbert
        float regularization_hilbert_strength_;
        /// whether to perform bidirectional moco
        bool bidirectional_moco_;

        /// image dissimilarity for cross-row moco
        Gadgetron::GT_IMAGE_DISSIMILARITY dissimilarity_cross_row_;
        /// number of resolution levels for cross-row moco
        unsigned int level_cross_row_;
        /// iterations for every level for cross-row moco
        std::vector<unsigned int> iter_cross_row_;
        /// regularization strength, hilbert for cross-row moco
        float regularization_hilbert_strength_cross_row_;
        /// whether to perform bidirectional moco for cross-row moco
        bool bidirectional_moco_cross_row_;
        /// whether to apply divergence free constraint for cross-row moco
        bool divergence_free_constraint_cross_row_;

        /// threshold for dissimilarity
        float dissimilarity_thres_;
        /// number of subdivision
        unsigned int div_num_;
        /// number of bidirectional iterations
        unsigned int inverse_deform_enforce_iter_;
        /// weight for bidirectional iteration
        float inverse_deform_enforce_weight_;

        /// percentage of images kept for average
        double percentage_kept_for_averaging_;

        /// whether to use soft averaging
        bool soft_averaging_;

        /// segmentation mode
        long segmentation_mode_;

        /// moco+averaging
        Gadgetron::GtMoCoAveraging<float, double, 2> mocoer2D_;

        ImageContainer2DMagType mag_container2D_;

        /// verbose mode for MOCO
        bool verboseModeMOCO_;
    };
}
