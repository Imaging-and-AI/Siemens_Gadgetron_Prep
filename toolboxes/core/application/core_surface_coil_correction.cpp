/** \file   core_surface_coil_correction.cpp
    \brief  Implement surface coil correction algorithms
    \author Hui Xue
*/

#include "core_surface_coil_correction.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "morphology.h"

namespace Gadgetron { 

    template<typename ValueType>
    GtSurfaceCoilCorrection<ValueType>::GtSurfaceCoilCorrection() : verbose_(false), performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);
    }

    template<typename ValueType>
    GtSurfaceCoilCorrection<ValueType>::
        ~GtSurfaceCoilCorrection()
    {
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image2DType& input, const Image2DType& gmap, size_t numOfRefinement, size_t gridSize[2], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image2DType& scc, hoNDArray<float>& mask)
    {
        try
        {
            scc = input;

            value_type totalResidual;

            FFD2DType ffd(input, gridSize[0], gridSize[1]);
            BSplineFFD<ValueType, float, 2, 1>& ffdBase = ffd;
            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Impl<2>(ffd, input, gmap, numOfRefinement, gridSize, useMask, noisebackground, thresRatioForNoise, scc, mask));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image2DType& input, ...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image3DType& input, const Image3DType& gmap, size_t numOfRefinement, size_t gridSize[3], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image3DType& scc, hoNDArray<float>& mask)
    {
        try
        {
            scc = input;

            value_type totalResidual;

            FFD3DType ffd(input, gridSize[0], gridSize[1], gridSize[2]);
            BSplineFFD<ValueType, float, 3, 1>& ffdBase = ffd;
            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Impl<3>(ffdBase, input, gmap, numOfRefinement, gridSize, useMask, noisebackground, thresRatioForNoise, scc, mask));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image3DType& input, ...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image4DType& input, const Image4DType& gmap, size_t numOfRefinement, size_t gridSize[4], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image4DType& scc, hoNDArray<float>& mask)
    {
        try
        {
            scc = input;

            FFD4DType ffd(input, gridSize[0], gridSize[1], gridSize[2], gridSize[3]);
            BSplineFFD<ValueType, float, 4, 1>& ffdBase = ffd;
            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Impl<4>(ffdBase, input, gmap, numOfRefinement, gridSize, useMask, noisebackground, thresRatioForNoise, scc, mask));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image4DType& input, ...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image2DType& input, size_t numOfRefinement, size_t gridSize[2], const MaskArrayType& mask, Image2DType& scc)
    {
        try
        {
            scc = input;

            FFD2DType ffd(input, gridSize[0], gridSize[1]);
            BSplineFFD<ValueType, float, 2, 1>& ffdBase = ffd;
            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Mask_Impl<2>(ffd, input, numOfRefinement, gridSize, mask, scc));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image2DType& input, size_t numOfRefinement, size_t gridSize[2], const MaskArrayType& mask, Image2DType& scc) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image3DType& input, size_t numOfRefinement, size_t gridSize[3], const MaskArrayType& mask, Image3DType& scc)
    {
        try
        {
            scc = input;

            FFD3DType ffd(input, gridSize[0], gridSize[1], gridSize[2]);
            BSplineFFD<ValueType, float, 3, 1>& ffdBase = ffd;
            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Mask_Impl<3>(ffdBase, input, numOfRefinement, gridSize, mask, scc));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image3DType& input, size_t numOfRefinement, size_t gridSize[3], const MaskArrayType& mask, Image3DType& scc) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image4DType& input, size_t numOfRefinement, size_t gridSize[4], const MaskArrayType& mask, Image4DType& scc)
    {
        try
        {
            scc = input;

            FFD4DType ffd(input, gridSize[0], gridSize[1], gridSize[2], gridSize[3]);
            BSplineFFD<ValueType, float, 4, 1>& ffdBase = ffd;
            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Mask_Impl<4>(ffdBase, input, numOfRefinement, gridSize, mask, scc));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD(const Image4DType& input, size_t numOfRefinement, size_t gridSize[4], const MaskArrayType& mask, Image4DType& scc) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(const Image2DType& input, const Image2DType& gmap, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[2], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image2DType& scc, hoNDArray<float>& mask)
    {
        try
        {
            typedef BSplineFFD<ValueType, float, 2, 1> FFDType;

            std::vector<FFDType*> ffds(numOfRefinementMax - numOfRefinementMin + 1, NULL);

            size_t f;
            for (f = 0; f<ffds.size(); f++)
            {
                ffds[f] = new FFD2DType(input, gridSize[0], gridSize[1]);
                ffds[f]->performTiming_ = this->performTiming_;
            }

            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Multiple_Impl<2>(ffds, input, gmap, numOfRefinementMin, numOfRefinementMax, gridSize, useMask, noisebackground, thresRatioForNoise, scc, mask));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(Image2DType&) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(const Image3DType& input, const Image3DType& gmap, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[3], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image3DType& scc, hoNDArray<float>& mask)
    {
        try
        {
            typedef BSplineFFD<ValueType, float, 3, 1> FFDType;

            std::vector<FFDType*> ffds(numOfRefinementMax - numOfRefinementMin + 1, NULL);

            size_t f;
            for (f = 0; f<ffds.size(); f++)
            {
                ffds[f] = new FFD3DType(input, gridSize[0], gridSize[1], gridSize[2]);
                ffds[f]->performTiming_ = this->performTiming_;
            }

            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Multiple_Impl<3>(ffds, input, gmap, numOfRefinementMin, numOfRefinementMax, gridSize, useMask, noisebackground, thresRatioForNoise, scc, mask));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(Image3DType) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(const Image4DType& input, const Image4DType& gmap, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[4], bool useMask, value_type noisebackground, value_type thresRatioForNoise, Image4DType& scc, hoNDArray<float>& mask)
    {
        try
        {
            typedef BSplineFFD<ValueType, float, 4, 1> FFDType;

            std::vector<FFDType*> ffds(numOfRefinementMax - numOfRefinementMin + 1, NULL);

            size_t f;
            for (f = 0; f<ffds.size(); f++)
            {
                ffds[f] = new FFD4DType(input, gridSize[0], gridSize[1], gridSize[2], gridSize[3]);
                ffds[f]->performTiming_ = this->performTiming_;
            }

            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Multiple_Impl<4>(ffds, input, gmap, numOfRefinementMin, numOfRefinementMax, gridSize, useMask, noisebackground, thresRatioForNoise, scc, mask));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(Image4DType) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(const Image2DType& input, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[2], const MaskArrayType& mask, Image2DType& scc)
    {
        try
        {
            typedef BSplineFFD<ValueType, float, 2, 1> FFDType;

            std::vector<FFDType*> ffds(numOfRefinementMax - numOfRefinementMin + 1, NULL);

            size_t f;
            for (f = 0; f<ffds.size(); f++)
            {
                ffds[f] = new FFD2DType(input, gridSize[0], gridSize[1]);
                ffds[f]->performTiming_ = this->performTiming_;
            }

            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Multiple_Mask_Impl<2>(ffds, input, numOfRefinementMin, numOfRefinementMax, gridSize, mask, scc));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(Image2DType, mask) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(const Image3DType& input, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[3], const MaskArrayType& mask, Image3DType& scc)
    {
        try
        {
            typedef BSplineFFD<ValueType, float, 3, 1> FFDType;

            std::vector<FFDType*> ffds(numOfRefinementMax - numOfRefinementMin + 1, NULL);

            size_t f;
            for (f = 0; f<ffds.size(); f++)
            {
                ffds[f] = new FFD3DType(input, gridSize[0], gridSize[1], gridSize[2]);
                ffds[f]->performTiming_ = this->performTiming_;
            }

            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Multiple_Mask_Impl<3>(ffds, input, numOfRefinementMin, numOfRefinementMax, gridSize, mask, scc));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(Image3DType, mask) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(const Image4DType& input, size_t numOfRefinementMin, size_t numOfRefinementMax, size_t gridSize[4], const MaskArrayType& mask, Image4DType& scc)
    {
        try
        {
            typedef BSplineFFD<ValueType, float, 4, 1> FFDType;

            std::vector<FFDType*> ffds(numOfRefinementMax - numOfRefinementMin + 1, NULL);

            size_t f;
            for (f = 0; f<ffds.size(); f++)
            {
                ffds[f] = new FFD4DType(input, gridSize[0], gridSize[1], gridSize[2], gridSize[3]);
                ffds[f]->performTiming_ = this->performTiming_;
            }

            GADGET_CHECK_RETURN_FALSE(this->computeSCC_FFD_Multiple_Mask_Impl<4>(ffds, input, numOfRefinementMin, numOfRefinementMax, gridSize, mask, scc));
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_FFD_Multiple(Image4DType, mask) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_LeastSquare(size_t RO, size_t E1,
        const T* PD, size_t boxFilterSize, value_type noisebackground, value_type thresRatioForNoise, T* scc)
    {
        T* PDFiltered = NULL;
        T* PDSquare = NULL;
        T* PDSquareFiltered = NULL;

        try
        {
            if (PD == NULL) return false;
            if (scc == NULL) return false;

            if (boxFilterSize % 2 == 0)
            {
                boxFilterSize++;
            }

            long long halfWidth = boxFilterSize / 2;

            // store filtered PD, PDSquare and PDSquareFiltered
            PDFiltered = new T[RO*E1];
            if (PDFiltered == NULL) return false;

            PDSquare = new T[RO*E1];
            if (PDSquare == NULL) return false;

            PDSquareFiltered = new T[RO*E1];
            if (PDSquareFiltered == NULL) return false;

            long long ro, e1;

            // compute PDSquare
            for (e1 = 0; e1<E1; e1++)
            {
                for (ro = 0; ro<RO; ro++)
                {
                    size_t offset = ro + e1 * RO;
                    PDSquare[offset] = PD[offset] * PD[offset];
                }
            }

            // box filtering PD and PDSquare
            T v = T(1.0 / (boxFilterSize*boxFilterSize));

            for (e1 = 0; e1<E1; e1++)
            {
                for (ro = 0; ro<RO; ro++)
                {
                    long long offset = ro + e1 * RO;

                    long long r, e, dr, de;

                    T vPD(0), vPDSquare(0);
                    for (r = -halfWidth; r <= halfWidth; r++)
                    {
                        dr = ro + r;
                        if (dr < 0) dr += RO;
                        if (dr >= RO) dr -= RO;

                        for (e = -halfWidth; e <= halfWidth; e++)
                        {
                            de = e1 + e;
                            if (de < 0) de += E1;
                            if (de >= E1) de -= E1;

                            vPD += PD[dr + de * RO] * v;
                            vPDSquare += PDSquare[dr + de * RO] * v;
                        }
                    }

                    PDFiltered[offset] = vPD;
                    PDSquareFiltered[offset] = vPDSquare;
                }
            }

            // compute the scc map
            v = noisebackground * thresRatioForNoise;
            for (e1 = 0; e1<E1; e1++)
            {
                for (ro = 0; ro<RO; ro++)
                {
                    long long offset = ro + e1 * RO;

                    scc[offset] = PDFiltered[offset] / (PDSquareFiltered[offset] + v * v);
                }
            }

            delete[] PDFiltered;
            delete[] PDSquare;
            delete[] PDSquareFiltered;
        }
        catch (...)
        {
            if (PDFiltered != NULL) delete[] PDFiltered;
            if (PDSquare != NULL) delete[] PDSquare;
            if (PDSquareFiltered != NULL) delete[] PDSquareFiltered;

            return false;
        }

        return true;
    }

    template<typename ValueType>
    bool GtSurfaceCoilCorrection<ValueType>::computeSCC_LeastSquare(const Image2DType& mag, const Image2DType& PD, const Image2DType& gmap, bool useMask, size_t boxFilterSize, value_type noisebackground, value_type thresRatioForNoise, Image2DType& scc, hoNDArray<float>& mask)
    {
        try
        {
            /*if (boxFilterSize % 2 == 0)
            {
            boxFilterSize++;
            }

            size_t RO = PD.get_size(0);
            size_t E1 = PD.get_size(1);

            scc = mag;

            if (!debugFolder_.empty()) gt_exporter_.export_image(PD, debugFolder_ + "PD");

            this->computeSCC_LeastSquare(RO, E1,
            PD.begin(), boxFilterSize, noisebackground, thresRatioForNoise, scc.begin());*/

            if (!debugFolder_.empty()) gt_exporter_.export_image(PD, debugFolder_ + "PD");

            hoNDArray<ValueType> ker(boxFilterSize, boxFilterSize);
            Gadgetron::fill(ker, ValueType(1));
            Gadgetron::scal((ValueType)(1.0 / (boxFilterSize*boxFilterSize)), ker);

            Image2DType PD_filtered(mag), PDPD_filtered(mag), PDPD;
            //Gadgetron::clear(PD_filtered);
            //Gadgetron::clear(PDPD_filtered);

            // Gadgetron::conv2(PD, ker, PD_filtered);
            // if (!debugFolder_.empty()) gt_exporter_.export_image(PD_filtered, debugFolder_ + "PD_filtered");

            Gadgetron::multiply(PD, PD, PDPD);
            if (!debugFolder_.empty()) gt_exporter_.export_image(PDPD, debugFolder_ + "PDPD");

            // Gadgetron::conv2(PDPD, ker, PDPD_filtered);
            // Gadgetron::addEpsilon(PDPD_filtered);
            // if (!debugFolder_.empty()) gt_exporter_.export_image(PDPD_filtered, debugFolder_ + "PDPD_filtered");

            if (boxFilterSize % 2 == 0)
            {
                boxFilterSize++;
            }
            long long halfWidth = boxFilterSize / 2;

            size_t RO = mag.get_size(0);
            size_t E1 = mag.get_size(1);

            T v = T(1.0 / (boxFilterSize*boxFilterSize));

            const T* pPD = PD.begin();
            T* pPDPD = PDPD.begin();
            T* pPD_filtered = PD_filtered.begin();
            T* pPDPD_filtered = PDPD_filtered.begin();

            long long ro, e1;
            for (e1 = 0; e1<E1; e1++)
            {
                for (ro = 0; ro<RO; ro++)
                {
                    long long offset = ro + e1 * RO;

                    long long r, e, dr, de;

                    T vPD(0), vPDSquare(0);
                    for (r = -halfWidth; r <= halfWidth; r++)
                    {
                        dr = ro + r;
                        if (dr < 0) dr = 0;
                        if (dr >= RO) dr = RO - 1;

                        for (e = -halfWidth; e <= halfWidth; e++)
                        {
                            de = e1 + e;
                            if (de < 0) de = 0;
                            if (de >= E1) de = E1 - 1;

                            vPD += pPD[dr + de * RO] * v;
                            vPDSquare += pPDPD[dr + de * RO] * v;
                        }
                    }

                    pPD_filtered[offset] = vPD;
                    pPDPD_filtered[offset] = vPDSquare;
                }
            }

            Gadgetron::addEpsilon(PDPD_filtered);

            if (!debugFolder_.empty()) gt_exporter_.export_image(PD_filtered, debugFolder_ + "PD_filtered");
            if (!debugFolder_.empty()) gt_exporter_.export_image(PDPD_filtered, debugFolder_ + "PDPD_filtered");

            scc = mag;

            std::vector<size_t> dim;
            mag.get_dimensions(dim);

            mask.create(dim);
            Gadgetron::clear(mask);

            value_type threshold = thresRatioForNoise * noisebackground;
            size_t N = scc.get_number_of_elements();

            bool hasGmap = true;
            size_t n;
            if (gmap.dimensions_equal(scc))
            {
                for (n = 0; n<N; n++)
                {
                    if (PD(n) > threshold*gmap(n))
                    {
                        mask(n) = 1;
                    }
                }
            }
            else
            {
                for (n = 0; n<N; n++)
                {
                    if (PD(n) > threshold)
                    {
                        mask(n) = 1;
                    }
                }

                hasGmap = false;
            }

            if (!debugFolder_.empty()) gt_exporter_.export_array(mask, debugFolder_ + "mask");

            hoNDArray<float> mask_hole_filling_filtered(mask);

            if (useMask)
            {
                // clean up the mask
                hoNDArray<unsigned int> label;
                bool isEightConnected = false;

                // ----------------------------------------------
                // hole filling for background

                label.clear();
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::bwlabel_2d(mask, float(0), label, isEightConnected));

                if (!debugFolder_.empty()) gt_exporter_.export_array(label, debugFolder_ + "mask_label");

                std::vector<unsigned int> labels, areas;
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::bwlabel_area_2d(label, labels, areas));

                size_t numLabels = labels.size();

                for (n = 0; n < N; n++)
                {
                    if (label(n)>0)
                    {
                        size_t ii;
                        for (ii = 0; ii < numLabels; ii++)
                        {
                            if (labels[ii] == label(n)) break;
                        }

                        if (ii < numLabels && areas[ii] < RO*E1 / 80.0)
                        {
                            mask(n) = 1;
                        }
                    }
                }
                if (!debugFolder_.empty()) gt_exporter_.export_array(mask, debugFolder_ + "mask_hole_filling");

                // hole filling for foreground

                isEightConnected = false;

                label.clear();
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::bwlabel_2d(mask, float(1), label, isEightConnected));

                if (!debugFolder_.empty()) gt_exporter_.export_array(label, debugFolder_ + "mask_label_foreground");

                labels.clear();
                areas.clear();
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::bwlabel_area_2d(label, labels, areas));

                numLabels = labels.size();

                for (n = 0; n < N; n++)
                {
                    if (label(n)>0)
                    {
                        size_t ii;
                        for (ii = 0; ii < numLabels; ii++)
                        {
                            if (labels[ii] == label(n)) break;
                        }

                        if (ii < numLabels && areas[ii] < RO*E1 / 320.0)
                        {
                            mask(n) = 0;
                        }
                    }
                }
                if (!debugFolder_.empty()) gt_exporter_.export_array(mask, debugFolder_ + "mask_hole_filling_foreground");

                // ----------------------------------------------

                // filtering the mask
                mask_hole_filling_filtered.copyFrom(mask);
                float sigma[2] = { 5.0f, 5.0f };
                Gadgetron::filterGaussian(mask_hole_filling_filtered, sigma);
                if (!debugFolder_.empty()) gt_exporter_.export_array(mask_hole_filling_filtered, debugFolder_ + "mask_hole_filling_filtered");
            }

            for (n = 0; n < N; n++)
            {
                double v1 = PD_filtered(n);
                double v2 = PDPD_filtered(n);

                if (useMask)
                {
                    if (mask(n) <= 0.5)
                    {
                        value_type v = thresRatioForNoise * noisebackground;
                        if (hasGmap) v *= gmap(n);
                        scc(n) = v1 / (v2 + v * v*(1.0 - mask_hole_filling_filtered(n)));
                    }
                    else
                    {
                        scc(n) = v1 / v2;
                    }
                }
                else
                {
                    value_type v = thresRatioForNoise * noisebackground;
                    if (hasGmap) v *= gmap(n);
                    scc(n) = v1 / (v2 + v * v);
                }
            }

            // Gadgetron::addEpsilon(scc);
            if (!debugFolder_.empty()) gt_exporter_.export_image(scc, debugFolder_ + "scc");

            Gadgetron::reciprocal_inplace(&scc);
            if (!debugFolder_.empty()) gt_exporter_.export_image(scc, debugFolder_ + "scc_reciprocal");
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GtSurfaceCoilCorrection<ValueType>::computeSCC_LeastSquare(Image2DType, mask) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType>
    void GtSurfaceCoilCorrection<ValueType>::print(std::ostream& os) const
    {
        using namespace std;

        os << "-------------- GT Surface coil correction -------------" << endl;
        os << "Differenc surface coil correction algorithm is implemented in this class " << endl;
        os << "Every calling function is independent from each other " << endl;
        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;
        os << "-------------------------------------------------------" << endl;
    }

    // ------------------------------------------------------------
    // Instantiation
    // ------------------------------------------------------------

    template class EXPORTGTTOOLBOXCORE GtSurfaceCoilCorrection< float >;
    template class EXPORTGTTOOLBOXCORE GtSurfaceCoilCorrection< double >;
}
