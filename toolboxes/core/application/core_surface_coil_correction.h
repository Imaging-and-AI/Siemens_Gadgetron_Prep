/** \file   GtSurfaceCoilCorrection.h
    \brief  Implement surface coil correction algorithms
    \author Hui Xue
*/

#pragma once


#include <gadgetron/hoMRImage.h>
#include <gadgetron/ho2DArray.h>
#include <gadgetron/ho3DArray.h>
#include <gadgetron/ho4DArray.h>
#include <gadgetron/ho5DArray.h>
#include <gadgetron/ho6DArray.h>
#include <gadgetron/ho7DArray.h>

#include <gadgetron/hoMatrix.h>
#include <gadgetron/hoNDArray_linalg.h>
#include <gadgetron/hoNDFFT.h>
#include <gadgetron/hoNDKLT.h>
#include <gadgetron/hoNDArray_utils.h>
#include <gadgetron/hoNDArray_elemwise.h>
#include <gadgetron/hoNDImage_util.h>
#include <gadgetron/hoMRImage.h>
#include <gadgetron/hoNDArray_reductions.h>
#include <gadgetron/hoNDArray_linalg.h>
#include <gadgetron/BSplineFFD2D.h>
#include <gadgetron/BSplineFFD3D.h>
#include <gadgetron/BSplineFFD4D.h>
#include <gadgetron/ImageIOAnalyze.h>

namespace Gadgetron { 

    template<typename ValueType>
    class GtSurfaceCoilCorrection
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

        typedef BSplineFFD2D<value_type, double, 1> FFD2DType;
        typedef BSplineFFD3D<value_type, double, 1> FFD3DType;
        typedef BSplineFFD4D<value_type, double, 1> FFD4DType;

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
    };
}
