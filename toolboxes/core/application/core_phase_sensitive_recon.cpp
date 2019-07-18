/** \file   GtPhaseSensitiveRecon.cpp
    \brief  Implement phase sensitive reconstruction
            The input is a 2D image container
    \author Hui Xue
*/

#include "core_phase_sensitive_recon.h"
#include <limits>
#include "hoNDFFT.h"
#include "hoNDImage_util.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"
#include "hoNDKLT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "hoMRImage.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_linalg.h"
#include "mri_core_kspace_filter.h"
#include "mri_core_def.h"

namespace Gadgetron { 

    template<typename ValueType, unsigned int D> 
    GtPhaseSensitiveRecon<ValueType, D>::GtPhaseSensitiveRecon() : performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        GADGET_CHECK_THROW(this->setDefaultParameters());
    }

    template<typename ValueType, unsigned int D> 
    GtPhaseSensitiveRecon<ValueType, D>::
    ~GtPhaseSensitiveRecon()
    {
    }

    template<typename ValueType, unsigned int D> 
    bool GtPhaseSensitiveRecon<ValueType, D>::setDefaultParameters()
    {
        row_PD_ = 1;
        col_compute_PSIR_windowing_ = 0;
        perform_PSIR_ = true;
        apply_PD_filtering_ = true;
        perform_SCC_PSIR_ = true;
        perform_SCC_MagIR_ = true;
        scale_factor_after_SCC_ = 250;
        offset_after_SCC_ = 0;
        compute_PSIR_windowing_ = true;

        scc_strategy_ = "FFDM";

        num_of_refinement_FFD_ = 3;
        num_of_refinement_max_FFD_ = 9;

        noise_masking_ = true;
        thres_ratio_noise_masking_ = 1.5f;

        preserve_PD_for_scc_ = true;

        size_t d;
        for ( d=0; d<D; d++ )
        {
            filter_width_[d] = 5;
        }

        intensity_scale_factor_ = 1;

        window_center_ = -1;
        window_width_ = -1;

        windowing_high_end_percentile_ = 0.95;

        verbose_ = false;

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool GtPhaseSensitiveRecon<ValueType, D>::
    initialize(const ImageContinerType& input)
    {
        try
        {
            size_t row = input.rows();
            std::vector<size_t> cols = input.cols();

            // allocate results
            std::vector<size_t> dim(row-1, cols[0]);
            GADGET_CHECK_RETURN_FALSE(PSIR_.create(dim));
            GADGET_CHECK_RETURN_FALSE(magIR_.create(dim));
            GADGET_CHECK_RETURN_FALSE(magIR_without_scc_.create(dim));
            GADGET_CHECK_RETURN_FALSE(PSIR_without_scc_.create(dim));

            dim.clear();
            dim.resize(1);
            dim[0] = cols[0];
            GADGET_CHECK_RETURN_FALSE(magPD_.create(dim));

            input_.clear();
            input_ = input;
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtPhaseSensitiveRecon<ValueType, D>::initialize(const TargetContinerType& targetContainer) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool GtPhaseSensitiveRecon<ValueType, D>::computeSCCMap(const ImageMagType& input, const ImageMagType& mag, ImageMagType& sccMap)
    {
        try
        {
            scc_.debugFolder_ = this->debugFolder_;

            bool useMask = this->noise_masking_;

            float scaleFactor = 1;

            try
            {
                scaleFactor = (float)input.attrib_.as_double(GADGETRON_IMAGE_SCALE_RATIO, 0);
            }
            catch(...)
            {
                useMask = false;
            }

            if ( !debugFolder_.empty() ) gt_exporter_.export_image(input, debugFolder_+"Input_For_SCC");

            sccMap = input;

            ImageMagType gmap;
            if (gmap_.dimensions_equal(input))
            {
                Gadgetron::complex_to_real(gmap_, gmap);
            }

            hoNDArray<float> mask;
            GADGET_CHECK_RETURN_FALSE(this->scc_.computeSCC_LeastSquare(mag, input, gmap, useMask, filter_width_[0], scaleFactor, thres_ratio_noise_masking_, sccMap, mask));
            if (!debugFolder_.empty()) gt_exporter_.export_image(sccMap, debugFolder_ + "SCC_LeastSquare");
            if (!debugFolder_.empty()) gt_exporter_.export_array(mask, debugFolder_ + "SCC_LeastSquare_mask");
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtPhaseSensitiveRecon<ValueType, D>::computeSCCMap(const ImageMagType& input, const ImageMagType& mag, ImageMagType& sccMap) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool GtPhaseSensitiveRecon<ValueType, D>::applyPDKSpaceFilter(ImageType& PDImage)
    {
        try
        {
            size_t RO = PDImage.get_size(0);
            size_t E1 = PDImage.get_size(1);

            if ( filter_PD_.get_size(0)!=RO || filter_PD_.get_size(1)!=E1 )
            {
                hoNDArray<ValueType> fRO(RO), fE1(E1);

                Gadgetron::generate_symmetric_filter(RO, fRO, ISMRMRD_FILTER_HANNING, 1.0, RO / 2);
                Gadgetron::generate_symmetric_filter(E1, fE1, ISMRMRD_FILTER_HANNING, 1.0, E1 / 2);

                Gadgetron::compute_2d_filter(fRO, fE1, filter_PD_);
            }

            // apply kspace filter
            hoNDArray<T> data(RO, E1, PDImage.begin());
            hoNDArray<T> dataFiltered(RO, E1);
            hoNDArray<T> dataKSpace(RO, E1);

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(data, dataKSpace);

            Gadgetron::apply_kspace_filter_ROE1(dataKSpace, filter_PD_, dataFiltered);

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataFiltered, dataKSpace);

            memcpy(PDImage.begin(), dataKSpace.begin(), sizeof(ValueType)*RO*E1);
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtPhaseSensitiveRecon<ValueType, D>::applyPDKSpaceFilter(ImageType& PDImage) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool GtPhaseSensitiveRecon<ValueType, D>::performPSIRRecon(const ImageContinerType& input)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize(input));

            size_t row = input_.rows();
            std::vector<size_t> cols = input_.cols();

            size_t RO = input_(0, 0).get_size(0);
            size_t E1 = input_(0, 0).get_size(1);

            size_t r, c;

            ImageType PDMag(RO, E1), PDPhase(RO, E1);

            for ( c=0; c<cols[0]; c++ )
            {
                ImageType& PDImage = input_(row_PD_, c);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PDImage, debugFolder_+"PDImage");

                // do not change magPD with filtering
                ImageType magPD;
                magPD.copyFrom(PDImage);
                Gadgetron::abs(PDImage, magPD);
                if (!debugFolder_.empty()) gt_exporter_.export_image_complex(magPD, debugFolder_ + "magPD");

                magPD_(0, c).copyFrom(magPD);

                if ( apply_PD_filtering_ )
                {
                    // apply the hanning filter
                    this->applyPDKSpaceFilter(PDImage);
                    if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PDImage, debugFolder_+"PDImageKSpaceFiltered");
                }

                Gadgetron::abs(PDImage, PDMag);
                Gadgetron::addEpsilon(PDMag);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PDMag, debugFolder_+"PDMag");

                if ( perform_PSIR_ )
                {
                    Gadgetron::divide(PDImage, PDMag, PDPhase);
                    if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PDPhase, debugFolder_+"PDPhase");

                    size_t rInd = 0;
                    for ( r=0; r<row; r++ )
                    {
                        if ( r == row_PD_ ) continue;

                        ImageType& IRImage = input_(r, c);
                        if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(IRImage, debugFolder_+"IRImage");

                        ImageType& magIRImage = magIR_(rInd, c);
                        // magIRImage = PSIRImage;
                        magIRImage = IRImage;

                        Gadgetron::abs(IRImage, magIRImage);
                        if (!debugFolder_.empty()) gt_exporter_.export_image_complex(magIRImage, debugFolder_ + "magIRImage");

                        magIR_without_scc_(rInd, c).copyFrom(magIRImage);

                        // ----------------------------

                        ImageType& PSIRImage = PSIR_(rInd, c);

                        PSIRImage = IRImage;
                        if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PSIRImage, debugFolder_+"PSIRImage");

                        Gadgetron::multiplyConj(PSIRImage, PDPhase, PSIRImage);
                        if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PSIRImage, debugFolder_+"PSIRImage_BackgroundPhaseRemoved");

                        Gadgetron::complex_to_real(PSIRImage);
                        if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PSIRImage, debugFolder_+"PSIRImage_real");

                        PSIR_without_scc_(rInd, c).copyFrom(PSIRImage);

                        rInd++;
                    }
                }
                else
                {
                    size_t rInd = 0;
                    for ( r=0; r<row; r++ )
                    {
                        if ( r == row_PD_ ) continue;

                        ImageType& IRImage = input_(r, c);

                        ImageType& magIRImage = magIR_(rInd, c);
                        magIRImage = IRImage;

                        Gadgetron::abs(IRImage, magIRImage);
                        if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(magIRImage, debugFolder_+"magIRImage");

                        magIR_without_scc_(rInd, c).copyFrom(magIRImage);

                        rInd++;
                    }
                }
            }

            if ( !debugFolder_.empty() )
            {
                hoNDArray<T> out;

                for ( r=0; r<row; r++ )
                {
                    input_.to_NDArray(r, out);
                    std::ostringstream ostr;
                    ostr << "Input_row" << r;
                    if ( !debugFolder_.empty() ) gt_exporter_.export_array_complex(out, debugFolder_+ostr.str());
                }

                if ( perform_PSIR_ )
                {
                    for ( r=0; r<row-1; r++ )
                    {
                        PSIR_.to_NDArray(r, out);
                        std::ostringstream ostr;
                        ostr << "PSIR_row" << r;
                        if ( !debugFolder_.empty() ) gt_exporter_.export_array_complex(out, debugFolder_+ostr.str());
                    }
                }

                for ( r=0; r<row-1; r++ )
                {
                    magIR_.to_NDArray(r, out);
                    std::ostringstream ostr;
                    ostr << "MagIR_row" << r;
                    if ( !debugFolder_.empty() ) gt_exporter_.export_array_complex(out, debugFolder_+ostr.str());
                }

                for ( r=0; r<row-1; r++ )
                {
                    magIR_without_scc_.to_NDArray(r, out);
                    std::ostringstream ostr;
                    ostr << "MagIR_without_scc_row" << r;
                    if ( !debugFolder_.empty() ) gt_exporter_.export_array_complex(out, debugFolder_+ostr.str());
                }

                magPD_.to_NDArray(0, out);
                std::ostringstream ostr;
                ostr << "magPD";
                if ( !debugFolder_.empty() ) gt_exporter_.export_array_complex(out, debugFolder_+ostr.str());
            }

            if (perform_SCC_PSIR_ || perform_SCC_MagIR_)
            {
                ImageType& magPD = magPD_(0, col_compute_PSIR_windowing_);

                // filter the magPD using median
                ImageMagType magPDReal(magPD.get_dimensions());
                Gadgetron::complex_to_real(magPD, magPDReal);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image(magPDReal, debugFolder_+"magPDReal");

                ImageMagType magPDFiltered(magPDReal);
                Gadgetron::filterMedian(magPDReal, filter_width_, magPDFiltered);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image(magPDFiltered, debugFolder_+"magPDFiltered");

                ImageType magPDFilteredCx(magPDFiltered.get_dimensions());
                Gadgetron::real_to_complex(magPDFiltered, magPDFilteredCx);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(magPDFilteredCx, debugFolder_+"magPDFilteredCx");

                ImageMagType sccMap(magPDReal);
                magPDReal.attrib_ = input(0, col_compute_PSIR_windowing_).attrib_;

                ImageType& magIRImage = magIR_(0, col_compute_PSIR_windowing_);
                ImageMagType magIRReal;
                Gadgetron::complex_to_real(magIRImage, magIRReal);

                GADGET_CHECK_RETURN_FALSE(this->computeSCCMap(magPDReal, magIRReal, sccMap));
                if ( !debugFolder_.empty() ) gt_exporter_.export_image(sccMap, debugFolder_+"sccMap");

                ImageType sccMapCx(sccMap.get_dimensions());
                Gadgetron::real_to_complex(sccMap, sccMapCx);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(sccMapCx, debugFolder_+"sccMapCx");

                ImageType PSIRImageBeforeSCC = PSIR_(0, col_compute_PSIR_windowing_);

                for ( c=0; c<cols[0]; c++ )
                {
                    for ( r=0; r<row-1; r++ )
                    {
                        if (perform_SCC_PSIR_)
                        {
                            ImageType& PSIRImage = PSIR_(r, c);
                            Gadgetron::divide(PSIRImage, sccMapCx, PSIRImage);
                            Gadgetron::scal(scale_factor_after_SCC_, PSIRImage);

                            if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(PSIRImage, debugFolder_+"PSIRImageSCC");
                        }

                        if (perform_SCC_MagIR_)
                        {
                            ImageType& magIRImage = magIR_(r, c);
                            Gadgetron::divide(magIRImage, sccMapCx, magIRImage);
                            Gadgetron::scal(scale_factor_after_SCC_, magIRImage);
                            if (!debugFolder_.empty()) gt_exporter_.export_image_complex(magIRImage, debugFolder_ + "magIRImageSCC");
                        }
                    }
                }

                // compute the windowing
                if ( perform_PSIR_ && compute_PSIR_windowing_ )
                {
                    ImageType& PSIRImage = PSIR_(0, col_compute_PSIR_windowing_);
                    ImageType& magIRImage = magIR_(0, col_compute_PSIR_windowing_);

                    this->calculate_window_level(magPDFilteredCx, PSIRImage, PSIRImageBeforeSCC, window_center_, window_width_);

                    GDEBUG_STREAM("PSIR - find windowing setting: " << window_center_ << " - " << window_width_);
                }

                for (c = 0; c<cols[0]; c++)
                {
                    for (r = 0; r<row - 1; r++)
                    {
                        if (perform_SCC_PSIR_)
                        {
                            ImageType& PSIRImage = PSIR_(r, c);
                            PSIRImage += offset_after_SCC_;
                            window_center_ += offset_after_SCC_;
                            if (!debugFolder_.empty()) gt_exporter_.export_image_complex(PSIRImage, debugFolder_ + "PSIRImageSCC");
                        }
                    }
                }

                if ( !debugFolder_.empty() )
                {
                    hoNDArray<T> out;

                    if ( perform_PSIR_ )
                    {
                        for ( r=0; r<row-1; r++ )
                        {
                            PSIR_.to_NDArray(r, out);
                            std::ostringstream ostr;
                            ostr << "PSIR_Norm_row" << r;
                            if ( !debugFolder_.empty() ) gt_exporter_.export_array_complex(out, debugFolder_+ostr.str());
                        }
                    }

                    for ( r=0; r<row-1; r++ )
                    {
                        magIR_.to_NDArray(r, out);
                        std::ostringstream ostr;
                        ostr << "MagIR_Norm_row" << r;
                        if ( !debugFolder_.empty() ) gt_exporter_.export_array_complex(out, debugFolder_+ostr.str());
                    }

                    std::ostringstream ostr;
                    ostr << "magPD_used_for_Norm";
                    if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(magPD, debugFolder_+ostr.str());

                    std::ostringstream ostr2;
                    ostr2 << "magPD_used_for_Norm_filtered";
                    if ( !debugFolder_.empty() ) gt_exporter_.export_image_complex(magPDFilteredCx, debugFolder_+ostr2.str());
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtPhaseSensitiveRecon<ValueType, D>::performPSIRRecon(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool GtPhaseSensitiveRecon<ValueType, D>::calculate_window_level(ImageType& magPDFiltered, ImageType& PSIRImage, ImageType& PSIRImageBeforeSCC, float& window_center, float& window_width)
    {
        try
        {
            try
            {
                intensity_scale_factor_ = (float)PSIRImage.attrib_.as_double(GADGETRON_IMAGE_SCALE_RATIO, 0);
            }
            catch(...)
            {
                intensity_scale_factor_ = 8;
            }

            // since gfactor is not taken into account here, we need higher thresholding ratio
            real_value_type thres = intensity_scale_factor_ * 5;

            size_t N = PSIRImage.get_number_of_elements();

            std::vector<size_t> dim;
            PSIRImage.get_dimensions(dim);
            ImageMagType mask(dim);
            Gadgetron::clear(mask);

            // get the foreground
            long long n;
            size_t numOfPixelInMask = 0;
            for ( n=0; n<N; n++ )
            {
                if ( magPDFiltered(n).real() > thres )
                {
                    mask(n) = 1;
                    numOfPixelInMask++;
                }
            }

            if ( !debugFolder_.empty() ) gt_exporter_.export_image(mask, debugFolder_+"PSIRMask");

            if ( numOfPixelInMask == 0 ) return true;

            // get the median of foreground
            std::vector<real_value_type> valueInMask(numOfPixelInMask, 0);
            size_t ind(0);
            for ( n=0; n<N; n++ )
            {
                if ( mask(n) == 1 )
                {
                    valueInMask[ind++] = PSIRImage(n).real();
                }
            }

            std::sort(valueInMask.begin(), valueInMask.end());
            real_value_type medianValueInMask = valueInMask[numOfPixelInMask/2];
            real_value_type w_high = valueInMask[ (size_t)(numOfPixelInMask *windowing_high_end_percentile_) ];

            //// make a 1D histogram on the mask of all foreground
            //size_t numOfBins = 256;
            //hoNDArray<double> hist1D(numOfBins), cdf1D(numOfBins);
            //hoNDArray<double> histBins(numOfBins);

            //Gadgetron::clear(hist1D);
            //Gadgetron::clear(cdf1D);
            //Gadgetron::clear(histBins);

            //real_value_type minValue = *std::min_element(valueInMask.begin(), valueInMask.end());
            //real_value_type maxValue = *std::max_element(valueInMask.begin(), valueInMask.end());

            //real_value_type range_t = real_value_type(1.0) / (maxValue - minValue + std::numeric_limits<real_value_type>::epsilon());

            //for (n = 0; n<numOfPixelInMask; n++)
            //{
            //    size_t indT = static_cast<size_t>(range_t*(valueInMask[n] - minValue)*(numOfBins - 1) + 0.5);
            //    hist1D(indT)++;
            //}

            //for (n = 0; n<numOfBins; n++)
            //{
            //    hist1D(n) /= numOfPixelInMask;
            //}

            //for (n = 0; n<numOfBins; n++)
            //{
            //    cdf1D(n) = 0;

            //    size_t m;
            //    for (m = 0; m <= n; m++)
            //    {
            //        cdf1D(n) += hist1D(m);
            //    }
            //}

            //long long n_high;
            //for (n_high = 0; n_high<numOfBins; n_high++)
            //{
            //    if (cdf1D(n_high) > windowing_high_end_percentile_) break;
            //}

            //long long n_low;
            //for (n_low = 0; n_low<numOfBins; n_low++)
            //{
            //    if (cdf1D(n_low) > 0.05) break;
            //}

            //real_value_type w_low = minValue + n_low * (maxValue - minValue) / numOfBins;
            //real_value_type w_high = minValue + n_high * (maxValue - minValue) / numOfBins;

            // get the second median and normal level, should be the center of myocardium
            ImageMagType thresdhold(dim);
            Gadgetron::clear(thresdhold);

            for ( n=0; n<N; n++ )
            {
                if ( (PSIRImage(n).real()<medianValueInMask) && (mask(n) == 1) )
                {
                    thresdhold(n) = 1;
                }
            }

            if ( !debugFolder_.empty() ) gt_exporter_.export_image(thresdhold, debugFolder_+"thresdhold");

            real_value_type normal_level = 0;

            numOfPixelInMask = 0;
            for ( n=0; n<N; n++ )
            {
                if ( thresdhold(n) == 1 )
                {
                    normal_level += PSIRImage(n).real();
                    numOfPixelInMask++;
                }
            }
            normal_level /= numOfPixelInMask;

            std::vector<real_value_type> valueInMaskMyo;
            valueInMaskMyo.resize(numOfPixelInMask, 0);

            ind = 0;
            for ( n=0; n<N; n++ )
            {
                if ( thresdhold(n) == 1 )
                {
                    valueInMaskMyo[ind++] = PSIRImage(n).real();
                }
            }
            std::nth_element(valueInMaskMyo.begin(), valueInMaskMyo.begin() + valueInMaskMyo.size() / 2, valueInMaskMyo.end());
            medianValueInMask = valueInMaskMyo[numOfPixelInMask / 2];

            // get the windowing setting

            real_value_type range = 1.1 * w_high - normal_level;
            real_value_type min_level = w_high - range;

            window_center = min_level + range / 2;
            window_width = range;

            //window_center = 2.1*medianValueInMask;
            //window_width = 2*(medianValueInMask-w_low);

            //if ( 2*(w_high-medianValueInMask) > window_width )
            //{
            //    window_width = 2*(w_high-medianValueInMask);
            //}

            //w_low = window_center - window_width / 2;

            //window_width *= 2.0f;
            //w_high = w_low + window_width;

            //window_center = (w_low + w_high) /2;

            //real_value_type f1dnormreal_val97 = minValue + n_high * (maxValue - minValue)/numOfBins;

            //real_value_type range = 1.1 * f1dnormreal_val97 - normal_level;
            //real_value_type min_level = f1dnormreal_val97 - range;

            //window_center = min_level + range/2;
            //window_width = range;
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtPhaseSensitiveRecon<ValueType, D>::calculate_window_level(ImageType& magPDFiltered, ImageType& PSIRImage, ImageType& PSIRImageBeforeSCC, float& window_center, float& window_width) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    void GtPhaseSensitiveRecon<ValueType, D>::print(std::ostream& os) const
    {
        using namespace std;

        os << "-------------- GTPlus PSIR Reconstruction -------------" << endl;

        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;

        os << "Row as the PD images is : " << row_PD_ << std::endl;
        os << "Column used for PSIR windowing computation is : " << col_compute_PSIR_windowing_ << std::endl;
        os << "-------------------------------------------------------" << endl;
    }

    // ------------------------------------------------------------
    // Instantiation
    // ------------------------------------------------------------

    template class EXPORTGTTOOLBOXCORE GtPhaseSensitiveRecon< std::complex<float>, 2 >;
    template class EXPORTGTTOOLBOXCORE GtPhaseSensitiveRecon< std::complex<double>, 2 >;
}
