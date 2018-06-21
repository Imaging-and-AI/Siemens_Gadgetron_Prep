/** \file   GtMoCoAveraging.h
    \brief  Implement 2D and 3D moco averaging
            The input is a 2D image container
            for every row, the moco+ave will be performed
            the different row can have the same or different reference frame
    \author Hui Xue
*/

#pragma once

#include "gadgetron_siemens_toolbox_core_export.h"
#include "hoImageRegContainer2DRegistration.h"
#include "mri_core_def.h"
#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"

namespace Gadgetron 
{ 
    template<typename ValueType, typename CoordType, unsigned int D> 
    class EXPORTGTTOOLBOXCORE GtMoCoAveraging
    {
    public:

        typedef GtMoCoAveraging<ValueType, CoordType, D> Self;

        typedef Gadgetron::hoMRImage<ValueType, D> ImageType;
        typedef Gadgetron::hoMRImage<ValueType, 2> Image2DType;
        typedef Gadgetron::hoMRImage<ValueType, 3> Image3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        typedef Gadgetron::hoImageRegContainer2DRegistration<ImageType, ImageType, CoordType> RegContainer2DType;
        typedef typename RegContainer2DType::TransformationDeformationFieldType TransformationDeformationFieldType;
        typedef typename RegContainer2DType::DeformationFieldType DeformationFieldType;

        typedef Gadgetron::hoNDImageContainer2D<ImageType> ImageContinerType;

        typedef Gadgetron::hoNDImageContainer2D<DeformationFieldType> DeformationFieldContinerType;

        typedef std::vector< std::vector< std::pair<double, size_t> > > MoCoQualityType;

        GtMoCoAveraging();
        virtual ~GtMoCoAveraging();

        /// set the default parameters
        virtual bool setDefaultParameters();

        /// compute moco and moco quality, but not averaging
        virtual bool computeMoCoAveraging(const ImageContinerType& input)
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize(input));
            GADGET_CHECK_RETURN_FALSE(this->pickReference());
            GADGET_CHECK_RETURN_FALSE(this->performMoCo());
            GADGET_CHECK_RETURN_FALSE(this->computeMoCoQuality());
            GADGET_CHECK_RETURN_FALSE(this->performAveraging(register_.warped_container_, averaged_));

            return true;
        }

        /// apply computed moco averaging on an input container and compute averaged images
        template<typename ImageType2> 
        bool applyMoCoAveraging(const Gadgetron::hoNDImageContainer2D< ImageType2 >& input, 
                Gadgetron::hoNDImageContainer2D< ImageType2 >& warpped, 
                Gadgetron::hoNDImageContainer2D< ImageType2 >& averaged)
        {
            GADGET_CHECK_RETURN_FALSE(register_.warpContainer2D(input, input, register_.deformation_field_, warpped));
            GADGET_CHECK_RETURN_FALSE(this->performAveraging(warpped, averaged));
            return true;
        }

        /// print the class information
        virtual void print(std::ostream& os) const;

        /// perform averaging
        template<typename ImageType2> 
        bool performAveraging( Gadgetron::hoNDImageContainer2D< ImageType2 >& input, 
                Gadgetron::hoNDImageContainer2D< ImageType2 >& averaged )
        {
            typedef typename ImageType2::value_type ValueType2;

            try
            {
                size_t row = input.rows();
                std::vector<size_t> cols = input.cols();

                std::vector<size_t> dim(1);
                dim[0] = row;
                GADGET_CHECK_RETURN_FALSE(averaged.create(dim));

                size_t r;
                for ( r=0; r<row; r++ )
                {
                    long long col = (long long)cols[r];
                    long long usedImageForAve = (long long)(col * percentage_kept_for_averaging_);
                    GDEBUG_STREAM("Number of used images for averaging : " << usedImageForAve);

                    size_t ref = reference_[r];

                    averaged(0, r) = input(r, ref);
                    if ( usedImageForAve <= 1 )
                    {
                        continue;
                    }
                    Gadgetron::clear(averaged(0, r));

                    long long ii;
                    for ( ii=0; ii<usedImageForAve; ii++ )
                    {
                        Gadgetron::add(averaged(0, r), input(r, moco_quality_measure_[r][ii].second ), averaged(0, r));

                        // mark the images picked for averaging
                        input(r, moco_quality_measure_[r][ii].second ).attrib_.append(GADGETRON_IMAGECOMMENT, "P");
                    }

                    Gadgetron::scal( (typename realType<ValueType2>::Type)( 1.0/std::sqrt( (double)usedImageForAve ) ),  averaged(0, r));
                }
            }
            catch(...)
            {
                GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::performAveraging(...) ... ");
                return false;
            }

            return true;
        }

        /// perform the cross-series moco on the averaged_ container
        virtual bool performCrossRowMoCo(size_t refRow);

        /// perform the cross-series moco on the averaged_ container
        /// every row can be specified a reference row
        /// if the reference of a row is itself, then no cross row moco is performed
        virtual bool performCrossRowMoCo(const std::vector<unsigned int>& refRows);

        template<typename ImageType2> 
        bool applyCrossRowMoCo(const Gadgetron::hoNDImageContainer2D< ImageType2 >& input, 
                Gadgetron::hoNDImageContainer2D< ImageType2 >& warpped)
        {
            GADGET_CHECK_RETURN_FALSE(register_cross_row_.warpContainer2D(input, input, register_cross_row_.deformation_field_, warpped));
            return true;
        }

        // ----------------------------------
        // parameters
        // ----------------------------------

        /// input image container containing images for processing
        ImageContinerType input_;

        /// the image container register
        RegContainer2DType register_;

        /// register for cross-row moco
        RegContainer2DType register_cross_row_;

        /// register for cross-row reference selection moco
        RegContainer2DType register_cross_row_ref_selection_;

        /// whether to perform cross row registration
        bool cross_row_reg_;

        /// force all rows in a group have the same reference
        bool cross_row_reference_;

        /// method to pick reference for every row in a group
        /// "SSD", "Deformation", "CrossCorrelation"
        /// if "Deformation" is used, then cross_row_reference_ is ignored
        std::string row_ref_pick_strategy_;

        /// if "Deformation" is used, indicate which row in a group is used as the target row
        size_t target_row_index_;

        /// method to compute moco quality
        /// "SSD", "Deformation"
        std::string moco_quality_strategy_;

        /// whether to determine moco quality in the ref picking phase
        /// if true, and row_ref_pick_strategy_ == "SSD", the SSD values will be used as moco quality
        /// if true, and row_ref_pick_strategy_ == "Deformation", the deformation field will be used as moco quality
        bool moco_quality_determined_in_ref_picking_;

        /// averaging percentage
        /// e.g., 0.5 means half of the data is used for averaging
        double percentage_kept_for_averaging_;

        /// whether to perform soft averaging
        bool soft_averaging_;

        /// whether need the mocoed images from input
        /// sometimes, only the deformation filed is needed
        bool moco_input_needed_;

        /// if multiple rows are given, indicate which rows belong to the same group
        /// by default, all rows belong to one group
        /// the cross_row_reference_ and cross_row_reg_ all 
        std::vector<size_t> row_group_index_;

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

        // ----------------------------------
        // registration results
        // ----------------------------------

        /// mocoed images and deformation fields are stored in the register_

        /// averaged images, [numOfRow 1]
        ImageContinerType averaged_;

        /// refrence picked for every row
        std::vector<unsigned int> reference_;

        /// moco quality measurements for every row
        /// the moco quality of a frame is measured by a double value
        MoCoQualityType moco_quality_measure_;

        bool initialize(const ImageContinerType& input);

        /// pick the reference
        virtual bool pickReference();

        /// perform moco averaging
        virtual bool performMoCo();

        /// compute moco quality measure
        bool computeMoCoQuality();

        /// compute moco quality measure using SSD
        bool computeMoCoQualitySSD(const ImageContinerType& warpped, const std::vector<unsigned int>& ref, MoCoQualityType& mocoQ);

        /// compute moco quality measure using deformation field
        bool computeMoCoQualityDeformation(const DeformationFieldContinerType* deform, const std::vector<unsigned int>& ref, MoCoQualityType& mocoQ);

        /// find row reference using median minimal SSD method
        bool findRefUsingSSD(size_t r, size_t& ind);

        /// find row reference using cross-correlation method
        bool findRefUsingCC(size_t r, size_t& ind, size_t start, size_t end);

        /// find row reference using deformtion method
        bool findRefUsingDeformation(size_t targetRow, size_t targetRef, const std::vector<size_t>& sourceRows, std::vector<size_t>& sourceRowRefs);
    };
}
