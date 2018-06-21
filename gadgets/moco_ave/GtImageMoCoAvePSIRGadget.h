/** \file   GtImageMoCoAvePSIRGadget.h
    \brief  The Gt image moco and averaging phase sensitive reconstruction gadget, used after Gt reconstruction for image data
            The PSIR dimension is SET
    \author Hui Xue
*/

#pragma once

#include "GtImageMoCoAveGadget.h"
#include "core_phase_sensitive_recon.h"

namespace Gadgetron { 

class EXPORTGTGADGET GtImageMoCoAvePSIRGadget : public GtImageMoCoAveGadget
{
public:
    GADGET_DECLARE(GtImageMoCoAvePSIRGadget);

    typedef GtImageMoCoAveGadget BaseClass;

    typedef BaseClass::ValueType ValueType;
    typedef BaseClass::Image2DType ImageType;
    typedef BaseClass::Image2DBufferType ImageBufferType;
    typedef BaseClass::ImgArrayType ImgArrayType;

    GtImageMoCoAvePSIRGadget();
    ~GtImageMoCoAvePSIRGadget();

    virtual int close(unsigned long flags);

    /// which set is the PD images
    unsigned int PD_set_;

    /// whether to apply extra filtering on PD
    bool apply_PD_filtering_;

    /// filter width for PD
    size_t filter_width_;

    /// surface coil correction methods
    /// "Median", "FFD", "FFDM"
    std::string scc_strategy_;

    /// for "FFD", number of refinements
    size_t num_of_refinement_FFD_;

    /// for "FFDM", number of maximal refinements
    /// the refinement range is [num_of_refinement_FFD_ num_of_refinement_max_FFD_]
    size_t num_of_refinement_max_FFD_;

    /// whether to perserve PD pixel values for scc map
    bool preserve_PD_for_scc_;

    /// whether to perform noise masking
    bool noise_masking_;

    /// threshold ratio for noise masking
    float thres_ratio_noise_masking_;

    /// scale factor after surface coil correction
    float scale_factor_after_SCC_;

    /// offset added to the PSIR
    float offset_after_SCC_;

    /// windowing high end percentile
    float windowing_high_end_percentile_;

    /// inversion time
    std::vector<float> TI_;

    /// whether to perform SCC on PSIR and MagIR
    bool perform_scc_PSIR_;
    bool perform_scc_mag_IR_;

    /// ori images
    bool send_ori_PSIR_;
    bool send_ori_mag_IR_;
    bool send_ori_mag_PD_;

    /// moco images
    bool send_moco_PSIR_;
    bool send_moco_mag_IR_;
    bool send_moco_mag_PD_;

    /// moco ave images
    bool send_moco_ave_PSIR_;
    bool send_moco_ave_mag_IR_;
    bool send_moco_ave_mag_PD_;

    /// whether to send mag IR without scc
    bool send_no_scc_mag_IR_;

    /// whether to send psir IR without scc
    bool send_no_scc_PSIR_;

    /// PSIR debug folder
    std::string debugFolder_PSIR_;
    std::string debugFolder_PSIR_fullPath_;

protected:

    GADGET_PROPERTY(debugFolder_PSIR, std::string, "Debug folder for psir only ... ", "");

    GADGET_PROPERTY(PD_set, int, "Indicate the set index of PD images", 1);
    GADGET_PROPERTY(apply_PD_filtering, bool, "Whether to apply extra filtering on PD images", true);
    GADGET_PROPERTY(preserve_PD_for_scc, bool, "Whether to use PD pixel values in scc map", true);
    GADGET_PROPERTY(filter_width, int, "Filter width for PD images", 7);
    GADGET_PROPERTY(scale_factor_after_SCC, double, "Scaling factor after surface coil correction", 1000);
    GADGET_PROPERTY(offset_after_SCC, double, "Offset to add after surface coil correction", 4096);

    GADGET_PROPERTY(windowing_high_end_percentile, double, "PSIR auto windowing, high end for histogram", 0.95);

    GADGET_PROPERTY_LIMITS(scc_strategy, std::string, "Strategy for surface coil correction", "LeastSquare",
        GadgetPropertyLimitsEnumeration, "Median", "FFD", "FFDM", "LeastSquare");

    GADGET_PROPERTY(num_of_refinement_FFD, int, "for 'FFD', number of refinements", 2);
    GADGET_PROPERTY(num_of_refinement_max_FFD, int, "for 'FFDM', number of maximal refinements", 7);

    GADGET_PROPERTY(noise_masking, bool, "Whether to apply noise mask in scc map", false);
    GADGET_PROPERTY(thres_ratio_noise_masking, double, "Threshold for masking the noise background", 3);

    GADGET_PROPERTY(perform_scc_PSIR, bool, "Whether to perform surface coil correction on PSIR images", true);
    GADGET_PROPERTY(perform_scc_mag_IR, bool, "Whether to perform surface coil correction on MagIR images", true);

    // ------------------------------------------------------------

    GADGET_PROPERTY(send_ori_PSIR, bool, "Whether to send original PSIR images", false);
    GADGET_PROPERTY(send_ori_mag_IR, bool, "Whether to send original MagIR images", false);
    GADGET_PROPERTY(send_ori_mag_PD, bool, "Whether to send original PD images", false);

    // ------------------------------------------------------------

    GADGET_PROPERTY(send_moco_PSIR, bool, "Whether to send moco PSIR images", false);
    GADGET_PROPERTY(send_moco_mag_IR, bool, "Whether to send moco MagIR images", false);
    GADGET_PROPERTY(send_moco_mag_PD, bool, "Whether to send moco PD images", false);

    // ------------------------------------------------------------

    GADGET_PROPERTY(send_moco_ave_PSIR, bool, "Whether to send moco ave PSIR images", true);
    GADGET_PROPERTY(send_moco_ave_mag_IR, bool, "Whether to send moco ave MagIR images", true);
    GADGET_PROPERTY(send_moco_ave_mag_PD, bool, "Whether to send moco ave PD images", false);

    // ------------------------------------------------------------

    GADGET_PROPERTY(send_no_scc_mag_IR, bool, "Whether to send MagIR images without surface coil correction, if perform_scc_mag_IR==true", false);
    GADGET_PROPERTY(send_no_scc_PSIR, bool, "Whether to send PSIR images without surface coil correction", false);

    // ------------------------------------------------------------

    // read in parameters
    bool readParameters();

    virtual int process_config(ACE_Message_Block* mb);
    virtual int processImageBuffer(ImageBufferType& ori);

    /// perform the PSIR
    bool performPSIR(ImageBufferType& input, ImageBufferType& magIR, ImageBufferType& magIRNoSCC, ImageBufferType& psIR, ImageBufferType& psIRNoSCC, ImageBufferType& magPD, std::vector<float>& wCenter, std::vector<float>& wWidth);
    bool performPSIR(ImageContainer2DType& input, ImageType& gmap, ImageContainer2DType& magIR, ImageContainer2DType& magIRNoSCC, ImageContainer2DType& psIR, ImageContainer2DType& psIRNoSCC, ImageContainer2DType& magPD, float& windowCenter, float& windowWidth);

    Gadgetron::GtPhaseSensitiveRecon<ValueType, 2> psir_reconer_;
};

}
