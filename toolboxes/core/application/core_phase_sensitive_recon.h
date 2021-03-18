/** \file   GtPhaseSensitiveRecon.h
    \brief  Implement phase sensitive reconstruction
            The input is a 2D image container
    \author Hui Xue
*/

#pragma once

#include "core_surface_coil_correction.h"
#include <ho2DArray.h>
#include <ho3DArray.h>
#include <ho4DArray.h>
#include <ho5DArray.h>
#include <ho6DArray.h>
#include <ho7DArray.h>
#include <hoMatrix.h>
#include <hoNDArray_linalg.h>
#include <hoNDFFT.h>
#include <hoNDKLT.h>
#include <hoNDArray_utils.h>
#include <hoNDArray_elemwise.h>
#include <hoNDImage_util.h>
#include <hoMRImage.h>
#include <hoNDArray_reductions.h>
#include <hoNDArray_linalg.h>
#include <hoNDImageContainer2D.h>

namespace Gadgetron { 

    template<typename ValueType, unsigned int D> 
    class GtPhaseSensitiveRecon
    {
    public:

        typedef GtPhaseSensitiveRecon<ValueType, D> Self;

        typedef Gadgetron::hoMRImage<ValueType, D> ImageType;
        typedef Gadgetron::hoMRImage<ValueType, 2> Image2DType;
        typedef Gadgetron::hoMRImage<ValueType, 3> Image3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef typename realType<T>::Type real_value_type;
        typedef Gadgetron::hoMRImage<real_value_type, D> ImageMagType;

        typedef Gadgetron::hoNDImageContainer2D<ImageType> ImageContinerType;

        typedef GtSurfaceCoilCorrection<real_value_type> SCCType;

        GtPhaseSensitiveRecon();
        virtual ~GtPhaseSensitiveRecon();

        /// set the default parameters
        virtual bool setDefaultParameters();

        /// compute the phase sensitive recon
        virtual bool performPSIRRecon(const ImageContinerType& input);

        /// compute the surface coil correction map
        bool computeSCCMap(const ImageMagType& input, const ImageMagType& mag, ImageMagType& sccMap);

        /// print the class information
        virtual void print(std::ostream& os) const;

        // ----------------------------------
        // parameters
        // ----------------------------------

        /// input image container containing images for processing
        ImageContinerType input_;

        /// gfactor map
        ImageType gmap_;

        /// row as the PD images
        unsigned int row_PD_;

        /// column used for PSIR windowing computation
        unsigned int col_compute_PSIR_windowing_;

        /// whether to perform PSIR
        bool perform_PSIR_;

        /// whether to apply filtering on the PD images before PSIR recon
        /// if true, an extra hanning filter is applied
        bool apply_PD_filtering_;

        /// whether to perform surface coil correction using PD
        bool perform_SCC_PSIR_;

        /// whether to perform SCC on magIR images
        bool perform_SCC_MagIR_;

        /// surface coil correction methods
        /// "Median", "FFD", "FFDM", "LeastSquare"
        std::string scc_strategy_;

        /// for "Median", median filter window width
        size_t filter_width_[D];

        /// for "FFD", number of refinements
        size_t num_of_refinement_FFD_;

        /// for "FFDM", number of maximal refinements
        /// the refinement range is [num_of_refinement_FFD_ num_of_refinement_max_FFD_]
        size_t num_of_refinement_max_FFD_;

        /// whether to perform noise masking
        bool noise_masking_;

        /// threshold ratio for noise masking
        float thres_ratio_noise_masking_;

        /// whether to preserve PD pixel values used for scc
        /// if false, the smoothed approximation will be used as scc map
        bool preserve_PD_for_scc_;

        /// scaling factor after SCC
        real_value_type scale_factor_after_SCC_;
        real_value_type offset_after_SCC_;

        /// whether to compute PSIR windowing
        bool compute_PSIR_windowing_;

        /// intensity scale factor for input images
        float intensity_scale_factor_;

        /// output parameters
        /// PSIR, magIR and magPD results
        ImageContinerType PSIR_;

        /// input image container containing images for processing
        ImageContinerType magIR_;

        /// in case of perform scc on magIR, this container has the non-scc magIR
        ImageContinerType magIR_without_scc_;

        /// this container has the non-scc PSIR
        ImageContinerType PSIR_without_scc_;

        /// input image container containing images for processing
        ImageContinerType magPD_;

        /// windowing 
        float windowing_high_end_percentile_;
        float window_center_;
        float window_width_;

        /// verbose mode
        bool verbose_;

        // ----------------------------------
        // debug and timing
        // ----------------------------------
        // clock for timing
        Gadgetron::GadgetronTimer gt_timer1_;
        Gadgetron::GadgetronTimer gt_timer2_;
        Gadgetron::GadgetronTimer gt_timer3_;

        bool performTiming_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // debug folder
        std::string debugFolder_;

    protected:

        /// surface coil correction
        SCCType scc_;

        /// store the kspace filter for PD images
        hoNDArray<ValueType> filter_PD_;

        bool initialize(const ImageContinerType& input);

        /// compute a proper window level/width for PSIR images
        bool calculate_window_level(ImageType& magPDFiltered, ImageType& PSIRImage, ImageType& PSIRImageBeforeSCC, float& window_center, float& window_width);

        /// apply PD kspace filter
        bool applyPDKSpaceFilter(ImageType& PDImage);
    };
}
