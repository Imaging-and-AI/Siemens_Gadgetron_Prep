/** \file   GtSurfaceCoilCorrection.h
    \brief  Implement surface coil correction algorithms
    \author Hui Xue
*/

#pragma once

#include "gadgetron_siemens_toolbox_core_export.h"
#include "hoMRImage.h"
#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"
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
#include "BSplineFFD2D.h"
#include "BSplineFFD3D.h"
#include "BSplineFFD4D.h"
#include "ImageIOAnalyze.h"

namespace Gadgetron { 

    template<typename ValueType>
    class EXPORTGTTOOLBOXCORE GtSurfaceCoilCorrection
    {
    public:

        typedef GtSurfaceCoilCorrection<ValueType> Self;

        typedef Gadgetron::hoMRImage<ValueType, 2> Image2DType;
        typedef Gadgetron::hoMRImage<ValueType, 3> Image3DType;
        typedef Gadgetron::hoMRImage<ValueType, 4> Image4DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;
        typedef typename realType<T>::Type real_value_type;

        typedef BSplineFFD2D<value_type, float, 1> FFD2DType;
        typedef BSplineFFD3D<value_type, float, 1> FFD3DType;
        typedef BSplineFFD4D<value_type, float, 1> FFD4DType;

        typedef hoNDArray<float> MaskArrayType;

        GtSurfaceCoilCorrection();
        virtual ~GtSurfaceCoilCorrection();

        /// compute the surface coil inhomogeneity map from the input image
        /// using the simple median filter approach
        /// kx, ky: the median filter kernel size
        template<unsigned int D>
        bool computeSCC_MedianFilter(const hoMRImage<ValueType, D>& input, size_t kSize[], hoMRImage<ValueType, D>& scc)
        {
            try
            {
                scc = input;
                Gadgetron::filterMedian(input, kSize, scc);
                Gadgetron::addEpsilon(scc);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image(scc, debugFolder_+"GtSurfaceCoilCorrection_MedianFilter_SCC");
            }
            catch(...)
            {
                GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_MedianFilter(hoMRImage<ValueType, D>) ... ");
                return false;
            }

            return true;
        }

        /// compute the surface coil map using the FFD approximations
        /// gridSize: the initial grid size for FFD
        /// if useMask is true, the noise background is firstly masked out
        /// noisebackground: the std of noise background, the threshold used to mask out the noise backgourd is thresRatioForNoise*noisebackground
        virtual bool computeSCC_FFD(const Image2DType& input, const Image2DType& gmap, size_t numOfRefinement, size_t gridSize[2], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image2DType& scc, hoNDArray<float>& mask);
        virtual bool computeSCC_FFD(const Image3DType& input, const Image3DType& gmap, size_t numOfRefinement, size_t gridSize[3], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image3DType& scc, hoNDArray<float>& mask);
        virtual bool computeSCC_FFD(const Image4DType& input, const Image4DType& gmap, size_t numOfRefinement, size_t gridSize[4], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image4DType& scc, hoNDArray<float>& mask);

        virtual bool computeSCC_FFD(const Image2DType& input, size_t numOfRefinement, size_t gridSize[2], const MaskArrayType& mask, Image2DType& scc);
        virtual bool computeSCC_FFD(const Image3DType& input, size_t numOfRefinement, size_t gridSize[3], const MaskArrayType& mask, Image3DType& scc);
        virtual bool computeSCC_FFD(const Image4DType& input, size_t numOfRefinement, size_t gridSize[4], const MaskArrayType& mask, Image4DType& scc);

        /// compute the surface coil map using the multiple FFD approximations
        /// sx, sy: the initial grid size for the FFD
        /// if useMask is true, the noise background is firstly masked out
        /// noisebackground: the std of noise background, the threshold used to mask out the noise backgourd is thresRatioForNoise*noisebackground
        /// numOfRefinementMin and numOfRefinementMax is the number of minimal and maximal refinement
        virtual bool computeSCC_FFD_Multiple(const Image2DType& input, const Image2DType& gmap, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[2], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image2DType& scc, hoNDArray<float>& mask);
        virtual bool computeSCC_FFD_Multiple(const Image3DType& input, const Image3DType& gmap, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[3], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image3DType& scc, hoNDArray<float>& mask);
        virtual bool computeSCC_FFD_Multiple(const Image4DType& input, const Image4DType& gmap, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[4], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image4DType& scc, hoNDArray<float>& mask);

        virtual bool computeSCC_FFD_Multiple(const Image2DType& input, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[2], const MaskArrayType& mask, Image2DType& scc);
        virtual bool computeSCC_FFD_Multiple(const Image3DType& input, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[3], const MaskArrayType& mask, Image3DType& scc);
        virtual bool computeSCC_FFD_Multiple(const Image4DType& input, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[4], const MaskArrayType& mask, Image4DType& scc);

        /// compute the surface coil map using the least square filter
        virtual bool computeSCC_LeastSquare(const Image2DType& mag, const Image2DType& PD, const Image2DType& gmap, bool useMask, size_t boxFilterSize, value_type noisebackground, value_type thresRatioForNoise, Image2DType& scc, hoNDArray<float>& mask);

        /// a simple version with only pure points
        virtual bool computeSCC_LeastSquare(size_t RO, size_t E1, const T* PD, size_t boxFilterSize, value_type noisebackground, value_type thresRatioForNoise, T* scc);

        /// print the class information
        virtual void print(std::ostream& os) const;

        // ----------------------------------
        // parameters
        // ----------------------------------

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

        template<unsigned int D>
        bool computeSCC_FFD_Mask_Impl(BSplineFFD<ValueType, float, D, 1>& ffd, const hoMRImage<ValueType, D>& input, size_t numOfRefinement, size_t gridSize[], const MaskArrayType& mask, hoMRImage<ValueType, D>& scc)
        {
            try
            {
                typedef hoMRImage<ValueType, D> ImageType;

                scc = input;

                value_type totalResidual;

                if ( !debugFolder_.empty() ) gt_exporter_.export_array(mask, debugFolder_+"GtSurfaceCoilCorrection_SCC_FFD_Mask");

                if ( this->performTiming_ ) { gt_timer1_.start("SCC FFD, approximation with mask ... "); }
                GADGET_CHECK_RETURN_FALSE(ffd.ffdApproxImage( const_cast<ImageType*>(&input), mask, totalResidual, numOfRefinement));
                if ( this->performTiming_ ) { gt_timer1_.stop(); }

                Gadgetron::addEpsilon(scc);
                if ( !debugFolder_.empty() ) gt_exporter_.export_image(scc, debugFolder_+"GtSurfaceCoilCorrection_SCC_FFD");
            }
            catch(...)
            {
                GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Mask_Impl(mask) ... ");
                return false;
            }

            return true;
        }

        template<unsigned int D>
        bool computeSCC_FFD_Impl(BSplineFFD<ValueType, float, D, 1>& ffd, const hoMRImage<ValueType, D>& input, const Image2DType& gmap, size_t numOfRefinement, size_t gridSize[], bool useMask, value_type noisebackground, value_type thresRatioForNoise, hoMRImage<ValueType, D>& scc, hoNDArray<float>& mask)
        {
            try
            {
                typedef hoMRImage<ValueType, D> ImageType;

                scc = input;

                if ( useMask )
                {
                    std::vector<size_t> dim;
                    input.get_dimensions(dim);

                    mask.create(dim);
                    Gadgetron::clear(mask);

                    value_type threshold = thresRatioForNoise*noisebackground;

                    size_t N = input.get_number_of_elements();

                    size_t n;
                    if (gmap.dimensions_equal(input))
                    {
                        for (n = 0; n<N; n++)
                        {
                            if (input(n) > threshold*gmap(n))
                            {
                                mask(n) = 1;
                            }
                        }
                    }
                    else
                    {
                        for (n = 0; n<N; n++)
                        {
                            if (input(n) > threshold)
                            {
                                mask(n) = 1;
                            }
                        }
                    }

                    if ( !debugFolder_.empty() ) gt_exporter_.export_array(mask, debugFolder_+"GtSurfaceCoilCorrection_SCC_FFD_Mask");

                    this->computeSCC_FFD_Mask_Impl(ffd, input, numOfRefinement, gridSize, mask, scc);
                }
                else
                {
                    value_type totalResidual;

                    ffd.performTiming_ = this->performTiming_;

                    if ( this->performTiming_ ) { gt_timer1_.start("SCC FFD, approximation ... "); }
                    GADGET_CHECK_RETURN_FALSE(ffd.ffdApproxImage( const_cast<ImageType&>(input), totalResidual, numOfRefinement));
                    if ( this->performTiming_ ) { gt_timer1_.stop(); }

                    if ( this->performTiming_ ) { gt_timer1_.start("SCC FFD, evaluate FFD ... "); }
                    GADGET_CHECK_RETURN_FALSE(ffd.evaluateFFDOnImage(scc));
                    if ( this->performTiming_ ) { gt_timer1_.stop(); }

                    Gadgetron::addEpsilon(scc);
                }

                if ( this->verbose_ )
                {
                    ffd.print(std::cout);
                }
            }
            catch(...)
            {
                GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Impl(...) ... ");
                return false;
            }

            return true;
        }

        template<unsigned int D>
        bool computeSCC_FFD_Multiple_Mask_Impl(std::vector< BSplineFFD<ValueType, float, D, 1>* >& ffds, const hoMRImage<ValueType, D>& input, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[], const MaskArrayType& mask, hoMRImage<ValueType, D>& scc)
        {
            try
            {
                typedef hoMRImage<ValueType, D> ImageType;

                scc = input;

                value_type totalResidual;
                size_t f;

                if ( this->performTiming_ ) { gt_timer1_.start("SCC FFD multiple, approximation with mask ... "); }
                for ( f=0; f<ffds.size(); f++ )
                {
                    GADGET_CHECK_RETURN_FALSE(ffds[f]->ffdApproxImage( const_cast<ImageType*>(&input), mask, totalResidual, numOfRefinementMin+f));
                }
                if ( this->performTiming_ ) { gt_timer1_.stop(); }

                Gadgetron::clear(scc);
                ImageType sccOneLevel(input);
                if ( this->performTiming_ ) { gt_timer1_.start("SCC FFD multiple, evaluate FFD multiple ... "); }
                for ( f=0; f<ffds.size(); f++ )
                {
                    GADGET_CHECK_RETURN_FALSE(ffds[f]->evaluateFFDOnImage(sccOneLevel));

                    std::ostringstream ostr;
                    ostr << "GtSurfaceCoilCorrection_SCC_FFD_Multiple_" << f;
                    if ( !debugFolder_.empty() ) gt_exporter_.export_image(sccOneLevel, debugFolder_+ostr.str());

                    Gadgetron::add(sccOneLevel, scc, scc);
                }
                if ( this->performTiming_ ) { gt_timer1_.stop(); }

                if ( this->verbose_ )
                {
                    for ( f=0; f<ffds.size(); f++ )
                    {
                        ffds[f]->print(std::cout);
                    }
                }

                for ( f=0; f<ffds.size(); f++ )
                {
                    delete ffds[f];
                }

                Gadgetron::addEpsilon(scc);
                Gadgetron::scal( (ValueType)(1.0/( ffds.size() )), scc );

                if ( !debugFolder_.empty() ) gt_exporter_.export_image(scc, debugFolder_+"GtSurfaceCoilCorrection_SCC_FFD_Multiple");
            }
            catch(...)
            {
                GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple_Mask_Impl(mask) ... ");
                return false;
            }

            return true;
        }

        template<unsigned int D>
        bool computeSCC_FFD_Multiple_Impl(std::vector< BSplineFFD<ValueType, float, D, 1>* >& ffds, const hoMRImage<ValueType, D>& input, const hoMRImage<ValueType, D>& gmap, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[], bool useMask, value_type noisebackground, value_type thresRatioForNoise, hoMRImage<ValueType, D>& scc, hoNDArray<float>& mask)
        {
            try
            {
                typedef hoMRImage<ValueType, D> ImageType;

                scc = input;

                value_type totalResidual;
                size_t f;

                if ( useMask )
                {
                    std::vector<size_t> dim;
                    input.get_dimensions(dim);
                    mask.create(dim);
                    Gadgetron::clear(mask);

                    value_type threshold = thresRatioForNoise*noisebackground;

                    size_t N = input.get_number_of_elements();

                    size_t n;

                    if (gmap.dimensions_equal(input))
                    {
                        for (n = 0; n<N; n++)
                        {
                            if (input(n) > threshold*gmap(n))
                            {
                                mask(n) = 1;
                            }
                        }
                    }
                    else
                    {
                        for (n = 0; n<N; n++)
                        {
                            if (input(n) > threshold)
                            {
                                mask(n) = 1;
                            }
                        }
                    }

                    if ( !debugFolder_.empty() ) gt_exporter_.export_array(mask, debugFolder_+"GtSurfaceCoilCorrection_SCC_FFD_Multiple_Mask");

                    GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Multiple_Mask_Impl(ffds, input, numOfRefinementMin, numOfRefinementMax, gridSize, mask, scc));
                }
                else
                {
                    if ( this->performTiming_ ) { gt_timer1_.start("SCC FFD multiple, approximation ... "); }
                    for ( f=0; f<ffds.size(); f++ )
                    {
                        GADGET_CHECK_RETURN_FALSE(ffds[f]->ffdApproxImage( const_cast<ImageType&>(input), totalResidual, numOfRefinementMin+f));
                    }
                    if ( this->performTiming_ ) { gt_timer1_.stop(); }

                    Gadgetron::clear(scc);
                    ImageType sccOneLevel(input);
                    if ( this->performTiming_ ) { gt_timer1_.start("SCC FFD multiple, evaluate FFD multiple ... "); }
                    for ( f=0; f<ffds.size(); f++ )
                    {
                        GADGET_CHECK_RETURN_FALSE(ffds[f]->evaluateFFDOnImage(sccOneLevel));

                        std::ostringstream ostr;
                        ostr << "GtSurfaceCoilCorrection_SCC_FFD_Multiple_" << f;
                        if ( !debugFolder_.empty() ) gt_exporter_.export_image(sccOneLevel, debugFolder_+ostr.str());

                        Gadgetron::add(sccOneLevel, scc, scc);
                    }
                    if ( this->performTiming_ ) { gt_timer1_.stop(); }

                    if ( this->verbose_ )
                    {
                        for ( f=0; f<ffds.size(); f++ )
                        {
                            ffds[f]->print(std::cout);
                        }
                    }

                    for ( f=0; f<ffds.size(); f++ )
                    {
                        delete ffds[f];
                    }

                    Gadgetron::addEpsilon(scc);
                    Gadgetron::scal( (ValueType)(1.0/( ffds.size() )), scc );
                    if ( !debugFolder_.empty() ) gt_exporter_.export_image(scc, debugFolder_+"GtSurfaceCoilCorrection_SCC_FFD_Multiple");
                    }
            }
            catch(...)
            {
                GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple_Impl(...) ... ");
                return false;
            }

            return true;
        }
    };
}
