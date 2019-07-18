/** \file   GtMoCoAveraging.cpp
    \brief  Implement 2D and 3D moco averaging
            The input is a 2D image container
            for every row, the moco+ave will be performed
            the different row can have the same or different reference frame
    \author Hui Xue
*/

#include "core_moco_averaging.h"
#include "hoImageRegContainer2DRegistration.h"

namespace Gadgetron { 

    struct compObj
    {
        compObj() {}
        ~compObj() {}

        bool operator()(const std::pair<double, size_t>& m1, const std::pair<double, size_t>& m2) const
        {
            return !(m1.first >= m2.first);
        }
    };

    template<typename ValueType, typename CoordType, unsigned int D> 
    GtMoCoAveraging<ValueType, CoordType, D>::GtMoCoAveraging() : performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        GADGET_CHECK_THROW(this->setDefaultParameters());
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    GtMoCoAveraging<ValueType, CoordType, D>::
    ~GtMoCoAveraging()
    {
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::setDefaultParameters()
    {
        GADGET_CHECK_RETURN_FALSE(register_.setDefaultParameters());

        cross_row_reg_ = false;
        cross_row_reference_ = false;
        row_ref_pick_strategy_ = "SSD";

        moco_quality_strategy_ = "Deformation";
        moco_quality_determined_in_ref_picking_ = true;

        percentage_kept_for_averaging_ = 0.5;
        soft_averaging_ = true;;
        moco_input_needed_ = true;

        target_row_index_ = 0;

        verbose_ = false;

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::
    initialize(const ImageContinerType& input)
    {
        try
        {
            size_t row = input.rows();
            std::vector<size_t> dim(1, row);

            // allocate averaged images
            GADGET_CHECK_RETURN_FALSE(averaged_.create(dim));

            reference_.resize(row, 0);

            moco_quality_measure_.resize(row);

            std::vector<size_t> cols = input.cols();

            size_t r, c;
            for ( r=0; r<row; r++ )
            {
                moco_quality_measure_[r].resize(cols[r]);
                for ( c=0; c<cols[r]; c++ )
                {
                    moco_quality_measure_[r][c].first = 0;
                    moco_quality_measure_[r][c].second = c;
                }
            }

            input_ = input;

            // all rows in one group
            if ( row_group_index_.size() != row ) row_group_index_.clear();
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::initialize(const TargetContinerType& targetContainer) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::findRefUsingSSD(size_t r, size_t& ind)
    {
        try
        {
            size_t row = input_.rows();
            std::vector<size_t> cols = input_.cols();

            GADGET_CHECK_RETURN_FALSE(r<row);

            long long col = (long long)cols[r];

            hoMatrixReal<ValueType> SSD(col, col);
            Gadgetron::clear(SSD);

            long long m, n;
            #pragma omp parallel default(none) private(m, n) shared(r, col, SSD)
            {
                ImageType diff(input_(r, 0).get_dimensions());

                #pragma omp for
                for ( m=0; m<col; m++ )
                {
                    ValueType v(0);
                    for ( n=m+1; n<col; n++ )
                    {
                        Gadgetron::subtract(input_(r, m), input_(r, n), diff);

                        v =Gadgetron::nrm2(diff);
                        SSD(m, n) = v;
                        SSD(n, m) = v;
                    }
                }

                diff.clear();
            }

            // sort for every column
            GADGET_CHECK_RETURN_FALSE(SSD.sort_ascending_along_row());

            // pick the middel row
            hoMatrixReal<ValueType> minimalSSDCol;
            GADGET_CHECK_RETURN_FALSE(SSD.subMatrix(minimalSSDCol, col/2, col/2, 0, col-1));

            // find minimal SSD
            ValueType minSSD(0);
            ValueType* pIm = minimalSSDCol.begin();

            size_t c;

            Gadgetron::minAbsolute(minimalSSDCol, minSSD, ind);

            /*minSSD = pIm[0];
            ind = 0;
            for (c=1; c<col; c++)
            {
                if (pIm[c]<=minSSD)
                {
                    minSSD= pIm[c];
                    ind = c;
                }
            }*/

            if ( ind >= (size_t)col ) ind=(size_t)col-1;

            if ( moco_quality_determined_in_ref_picking_ )
            {
                // fill in the moco_quality_measure_
                for ( c=0; c<col; c++ )
                {
                    moco_quality_measure_[r][c].first = (double)( minimalSSDCol(0, c) );
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in findRefUsingSSD(size_t r, size_t& ind) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::findRefUsingCC(size_t r, size_t& ind, size_t start, size_t end)
    {
        try
        {
            size_t row = input_.rows();
            std::vector<size_t> cols = input_.cols();

            GADGET_CHECK_RETURN_FALSE(r<row);

            long long col = (long long)cols[r];

            // compute std image
            ImageType stdImage = input_(r, 0);

            size_t N = input_(r, 0).get_number_of_elements();

            hoNDArray<T> buf(end-start+1);

            size_t ii, c;
            for ( ii=0; ii<N; ii++ )
            {
                for ( c=start; c<=end; c++ )
                {
                    buf[c-start] = input_(r, c)(ii);
                }
                stdImage(ii) = Gadgetron::stddev(&buf);
            }

            if ( !debugFolder_.empty() ) gt_exporter_.export_image(stdImage, debugFolder_+"stdImage");

            hoNDArray<T> CC(end-start+1);

            for ( c=start; c<=end; c++ )
            {
                typename ImageType::value_type cc;
                corrCoef(stdImage, input_(r, c), cc);
                CC(c-start) = cc;
            }

            T* pIm = CC.begin();
            T maxCC = pIm[0];
            ind = 0;
            for (c=1; c<end-start+1; c++)
            {
                if (pIm[c]>=maxCC)
                {
                    maxCC= pIm[c];
                    ind = c;
                }
            }

            ind += start;

            if ( ind >= (size_t)col ) ind=(size_t)col-1;

            if ( moco_quality_determined_in_ref_picking_ )
            {
                // fill in the moco_quality_measure_
                for ( c=start; c<=end; c++ )
                {
                    moco_quality_measure_[r][c].first = (double)( -CC(c-start) );
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in findRefUsingCC(size_t r, size_t& ind, size_t start, size_t end) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::findRefUsingDeformation(size_t targetRow, size_t targetRef, const std::vector<size_t>& sourceRows, std::vector<size_t>& sourceRowRefs)
    {
        try
        {
            size_t row = input_.rows();
            std::vector<size_t> cols = input_.cols();

            GADGET_CHECK_RETURN_FALSE(targetRow<row);
            GADGET_CHECK_RETURN_FALSE(targetRef<cols[targetRow]);

            size_t numOfSourceRows = sourceRows.size();
            sourceRowRefs.resize(numOfSourceRows, targetRef);

            bool warped = false;

            if ( this->performTiming_ ) { gt_timer3_.start("GtMoCoAveraging - perform cross-row moco for reference selection ... "); }

            size_t totalNumOfImages = 1;

            size_t r;
            for ( r=0; r<numOfSourceRows; r++ )
            {
                totalNumOfImages  += cols[ sourceRows[r] ];
            }

            std::vector<size_t> colOfMoCoArray(1, totalNumOfImages);

            ImageContinerType imageArrrayForMOCO;
            imageArrrayForMOCO.create(colOfMoCoArray, false);

            imageArrrayForMOCO.set( &input_(targetRow, targetRef), 0, 0);

            size_t ind(1);
            for ( r=0; r<numOfSourceRows; r++ )
            {
                size_t c;
                for ( c=0; c<cols[ sourceRows[r] ]; c++ )
                {
                    imageArrrayForMOCO.set( &input_(sourceRows[r], c), 0, ind++);
                }
            }

            std::vector<unsigned int> refMoCo(1, 0);

            if ( !debugFolder_.empty() )
            {
                warped = true;

                hoNDArray<T> out;
                imageArrrayForMOCO.to_NDArray(0, out);

                std::ostringstream ostr;
                ostr << "ImageArrrayForRefSelectionBeforeMOCO";

                if ( !debugFolder_.empty() ) gt_exporter_.export_array(out, debugFolder_+ostr.str());
            }

            register_cross_row_ref_selection_.registerOverContainer2DFixedReference(imageArrrayForMOCO, refMoCo, warped);

            if ( !debugFolder_.empty() )
            {
                warped = true;

                hoNDArray<T> out;
                register_cross_row_ref_selection_.warped_container_.to_NDArray(0, out);

                std::ostringstream ostr;
                ostr << "ImageArrrayForRefSelectionAfterMOCO";

                if ( !debugFolder_.empty() ) gt_exporter_.export_array(out, debugFolder_+ostr.str());
            }

            hoImageRegDeformationField<CoordType, D> deformTransform;
            CoordType meanDeform, maxDeform, meanLogJac, maxLogJac;
            hoNDArray<CoordType> jac;
            DeformationFieldType* deform_field[D];

            ind = 1;
            for ( r=0; r<numOfSourceRows; r++ )
            {
                size_t d, c;

                CoordType minMeanDeform = (CoordType)1e6;
                CoordType minMeanLogJac = (CoordType)1e6;
                size_t refSource=0;

                for ( c=0; c<cols[ sourceRows[r] ]; c++ )
                {
                    for ( d=0; d<D; d++ )
                    {
                        deform_field[d] = &register_cross_row_ref_selection_.deformation_field_[d](0, ind);
                    }
                    ind++;

                    deformTransform.jacobianPosition(jac, deform_field, 2);
                    deformTransform.analyzeJacobianAndDeformation(jac, deform_field, meanDeform, maxDeform, meanLogJac, maxLogJac, 2);

                    //if ( meanLogJac < minMeanLogJac )
                    //{
                    //    minMeanLogJac = meanLogJac;
                    //    refSource = c;
                    //}

                    if (meanDeform < minMeanDeform)
                    {
                        minMeanDeform = meanDeform;
                        refSource = c;
                    }

                    if ( moco_quality_determined_in_ref_picking_ )
                    {
                        // fill in the moco_quality_measure_
                        // moco_quality_measure_[ sourceRows[r] ][c].first = (double)meanLogJac;
                        moco_quality_measure_[sourceRows[r]][c].first = (double)meanDeform;
                    }
                }

                sourceRowRefs[r] = refSource;
            }

            if ( this->performTiming_ ) { gt_timer3_.stop(); }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in findRefUsingDeformation(size_t targetRow, size_t targetRef, const std::vector<size_t>& sourceRows, std::vector<size_t>& sourceRowRefs) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::pickReference()
    {
        try
        {
            size_t row = input_.rows();
            std::vector<size_t> cols = input_.cols();
            size_t r;
            size_t ind(0);

            if ( row_ref_pick_strategy_=="SSD" || row_ref_pick_strategy_=="CrossCorrelation" )
            {
                for ( r=0; r<row; r++ )
                {
                    if ( row_ref_pick_strategy_=="SSD" )
                    {
                        if ( !this->findRefUsingSSD(r, ind) )
                        {
                            ind = 0;
                        }
                    }
                    else
                    {
                        if ( !this->findRefUsingCC(r, ind, 0, cols[r]-1) )
                        {
                            ind = 0;
                        }
                    }

                    reference_[r] = (unsigned int)ind;
                }

                if ( cross_row_reference_ && row>1 )
                {
                    if ( row_group_index_.size() != row )
                    {
                        std::vector<unsigned int> refComputed(reference_);

                        std::sort(refComputed.begin(), refComputed.end());
                        size_t commonRef = refComputed[row/2];

                        for ( r=0; r<row; r++ )
                        {
                            if ( commonRef < cols[r] )
                            {
                                reference_[r] = (unsigned int)commonRef;
                            }
                        }
                    }
                    else
                    {
                        size_t g;
                        std::vector<size_t> groupInd;

                        for ( r=0; r<row; r++ )
                        {
                            bool newGroup = true;
                            for ( g=0; g<groupInd.size(); g++ )
                            {
                                if (groupInd[g]==row_group_index_[r])
                                {
                                    newGroup = false;
                                    break;
                                }
                            }

                            if ( newGroup ) groupInd.push_back(row_group_index_[r]);
                        }

                        for ( g=0; g<groupInd.size(); g++ )
                        {
                            std::vector<size_t> refInd;
                            for ( r=0; r<row; r++ )
                            {
                                if (groupInd[g]==row_group_index_[r])
                                {
                                    refInd.push_back(reference_[r]);
                                }
                            }

                            std::sort(refInd.begin(), refInd.end());
                            size_t commonRef = refInd[refInd.size()/2];

                            for ( r=0; r<row; r++ )
                            {
                                if (groupInd[g]==row_group_index_[r])
                                {
                                    reference_[r] = (unsigned int)commonRef;
                                }
                            }
                        }
                    }
                }
            }
            else if (row_ref_pick_strategy_ == "Deformation")
            {
                if ( row_group_index_.size() != row )
                {
                    // only one group is available

                    // find the reference for the target row
                    if ( !this->findRefUsingSSD(target_row_index_, ind) )
                    {
                        GERROR_STREAM("failed in findRefUsingSSD(target_row_index_, ind) ... ");
                        ind = 0;
                    }

                    reference_[target_row_index_] = (unsigned int)ind;

                    // images in other rows are registered to reference
                    std::vector<size_t> sourceRowInd, refOfSourceRows(row-1, 0);
                    for ( r=0; r<row; r++ )
                    {
                        if ( r == target_row_index_ ) continue;
                        sourceRowInd.push_back(r);
                    }

                    if ( !this->findRefUsingDeformation(target_row_index_, reference_[target_row_index_], sourceRowInd, refOfSourceRows) )
                    {
                        refOfSourceRows.resize(row-1, reference_[target_row_index_]);
                    }

                    ind=0;
                    for ( r=0; r<row; r++ )
                    {
                        if ( r == target_row_index_ ) continue;
                        reference_[r] = (unsigned int)refOfSourceRows[ind++];
                    }
                }
                else
                {
                    size_t g;
                    std::vector<size_t> groupInd;

                    for ( r=0; r<row; r++ )
                    {
                        bool newGroup = true;
                        for ( g=0; g<groupInd.size(); g++ )
                        {
                            if (groupInd[g]==row_group_index_[r])
                            {
                                newGroup = false;
                                break;
                            }
                        }

                        if ( newGroup ) groupInd.push_back(row_group_index_[r]);
                    }

                    for ( g=0; g<groupInd.size(); g++ )
                    {
                        std::vector<size_t> rowsInGroup;
                        for ( r=0; r<row; r++ )
                        {
                            if (groupInd[g]==row_group_index_[r])
                            {
                                rowsInGroup.push_back(r);
                            }
                        }

                        if ( !this->findRefUsingSSD(rowsInGroup[target_row_index_], ind) )
                        {
                            GERROR_STREAM("failed in findRefUsingSSD(rowsInGroup[target_row_index_], ind) ... ");
                            ind = 0;
                        }
                        reference_[ rowsInGroup[target_row_index_] ] = (unsigned int)ind;

                        std::vector<size_t> sourceRowInd, refOfSourceRows(rowsInGroup.size()-1, 0);
                        for ( r=0; r<rowsInGroup.size(); r++ )
                        {
                            if ( r == target_row_index_ ) continue;
                            sourceRowInd.push_back(rowsInGroup[r]);
                        }

                        if ( !this->findRefUsingDeformation(rowsInGroup[target_row_index_], reference_[ rowsInGroup[target_row_index_] ], sourceRowInd, refOfSourceRows) )
                        {
                            refOfSourceRows.resize(row-1, reference_[ rowsInGroup[target_row_index_] ]);
                        }

                        ind=0;
                        for ( r=0; r<rowsInGroup.size(); r++ )
                        {
                            if ( r == target_row_index_ ) continue;
                            reference_[ rowsInGroup[r] ] = (unsigned int)refOfSourceRows[ind++];
                        }
                    }
                }
            }

            if ( verbose_ )
            {
                GDEBUG_STREAM("GtMoCoAveraging - reference : ");

                if ( row_group_index_.size() != row )
                {
                    for ( r=0; r<row; r++ )
                    {
                        GDEBUG_STREAM("row " << r << " - " << reference_[r]);
                    }
                }
                else
                {
                    for ( r=0; r<row; r++ )
                    {
                        GDEBUG_STREAM("row " << r << " - group " << row_group_index_[r] << " - ref index " << reference_[r]);
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::pickReference(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::performMoCo()
    {
        try
        {
            register_.debugFolder_ = debugFolder_;
            register_.verbose_ = verbose_;

            if ( verbose_ )
            {
                register_.print(std::cout);
            }

            bool warped = moco_input_needed_;

            if ( this->performTiming_ ) { gt_timer3_.start("GtMoCoAveraging - perform moco ... "); }

            if ( register_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE )
            {
                register_.registerOverContainer2DFixedReference(input_, reference_, warped);
            }

            if ( register_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_PROGRESSIVE )
            {
                warped = true;
                register_.registerOverContainer2DProgressive(input_, reference_);
            }

            if ( this->performTiming_ ) { gt_timer3_.stop(); }

            if ( !debugFolder_.empty() && warped )
            {
                hoNDArray<T> out;

                size_t row = input_.rows();

                size_t r;
                for ( r=0; r<row; r++ )
                {
                    register_.warped_container_.to_NDArray(r, out);

                    std::ostringstream ostr;
                    ostr << "MOCO_row" << r;

                    if ( !debugFolder_.empty() ) gt_exporter_.export_array(out, debugFolder_+ostr.str());
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::performMoCo(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::performCrossRowMoCo(size_t refRow)
    {
        try
        {
            size_t col = averaged_.cols(0);
            GADGET_CHECK_RETURN_FALSE(refRow < col);

            register_cross_row_.debugFolder_ = debugFolder_;
            register_cross_row_.verbose_ = verbose_;

            if ( verbose_ )
            {
                register_cross_row_.print(std::cout);
            }

            bool warped = moco_input_needed_;

            if ( this->performTiming_ ) { gt_timer3_.start("GtMoCoAveraging - perform cross-row moco ... "); }

            if (register_cross_row_.apply_divergence_free_constraint_)
            {
                GDEBUG_STREAM("Cross row moco is with the divergence free constraint ... ");
            }

            if (!debugFolder_.empty() && warped)
            {
                hoNDArray<T> out;
                averaged_.to_NDArray(0, out);

                std::ostringstream ostr;
                ostr << "MOCO_cross_row_Input";

                if (!debugFolder_.empty()) gt_exporter_.export_array(out, debugFolder_ + ostr.str());
            }

            std::vector<unsigned int> refCrossRow(1, (unsigned int)refRow);

            if ( register_cross_row_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE )
            {
                register_cross_row_.registerOverContainer2DFixedReference(averaged_, refCrossRow, warped);
            }

            if ( register_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_PROGRESSIVE )
            {
                warped = true;
                register_cross_row_.registerOverContainer2DProgressive(averaged_, refCrossRow);
            }

            if ( this->performTiming_ ) { gt_timer3_.stop(); }

            if ( !debugFolder_.empty() && warped )
            {
                hoNDArray<T> out;
                register_cross_row_.warped_container_.to_NDArray(0, out);

                std::ostringstream ostr;
                ostr << "MOCO_cross_row";

                if ( !debugFolder_.empty() ) gt_exporter_.export_array(out, debugFolder_+ostr.str());
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::performCrossRowMoCo(size_t refRow) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::performCrossRowMoCo(const std::vector<unsigned int>& refRows)
    {
        try
        {
            size_t col = averaged_.cols(0);
            GADGET_CHECK_RETURN_FALSE(refRows.size() >= col);

            register_cross_row_.debugFolder_ = debugFolder_;
            register_cross_row_.verbose_ = verbose_;

            if ( verbose_ )
            {
                register_cross_row_.print(std::cout);
            }

            if (register_cross_row_.apply_divergence_free_constraint_)
            {
                GDEBUG_STREAM("Cross row moco is with the divergence free constraint ... ");
            }

            bool warped = moco_input_needed_;

            ImageContinerType targetAveraged;
            std::vector<size_t> colSize = averaged_.cols();
            targetAveraged.create(colSize, false);
            targetAveraged.delete_data_on_destruct(false);

            size_t ii;
            for ( ii=0; ii<col; ii++ )
            {
                targetAveraged.set( &averaged_(0, refRows[ii]), 0, ii );
            }

            if (!debugFolder_.empty() && warped)
            {
                hoNDArray<T> out;
                targetAveraged.to_NDArray(0, out);

                std::ostringstream ostr;
                ostr << "MOCO_cross_row_Input";

                if (!debugFolder_.empty()) gt_exporter_.export_array(out, debugFolder_ + ostr.str());
            }

            if ( this->performTiming_ ) { gt_timer3_.start("GtMoCoAveraging - perform cross-row moco ... "); }

            register_cross_row_.registerOverContainer2DPairWise(targetAveraged, averaged_, warped);

            if ( this->performTiming_ ) { gt_timer3_.stop(); }

            if ( !debugFolder_.empty() && warped )
            {
                hoNDArray<T> out;
                register_cross_row_.warped_container_.to_NDArray(0, out);

                std::ostringstream ostr;
                ostr << "MOCO_cross_row";

                if ( !debugFolder_.empty() ) gt_exporter_.export_array(out, debugFolder_+ostr.str());
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::performCrossRowMoCo(const std::vector<unsigned int>& refRows) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualitySSD(const ImageContinerType& warpped, const std::vector<unsigned int>& reference, MoCoQualityType& mocoQ)
    {
        try
        {
            size_t row = warpped.rows();
            std::vector<size_t> cols = warpped.cols();

            long long r;

            for ( r=0; r<(long long)row; r++ )
            {
                size_t ref = reference[r];
                mocoQ[r][ref].first = 0;

                long long c;

                ValueType sumSSD = 0;
                for ( c=0; c<(long long)cols[r]; c++ )
                {
                    sumSSD += mocoQ[r][c].first;
                }

                if ( !moco_quality_determined_in_ref_picking_ || std::abs(sumSSD)<0.001 )
                {
                    const ImageType& refImage = warpped(r, ref);

                    {
                        #pragma omp parallel default(none) private(c) shared(refImage, cols, r, warpped, mocoQ, ref) if (D==2)
                        {
                            ImageType diff(refImage);
                            ValueType v;

                            #pragma omp for 
                            for ( c=0; c<(long long)cols[r]; c++ )
                            {
                                if ( c == ref ) continue;

                                const ImageType& warppedImage = warpped(r, c);

                                Gadgetron::subtract(refImage, warppedImage, diff);
                                v = Gadgetron::nrm2(diff);
                                mocoQ[r][c].first = v;
                            }

                            diff.clear();
                        }
                    }
                }

                std::sort(mocoQ[r].begin(), mocoQ[r].end(), Gadgetron::compObj() );
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualitySSD(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualityDeformation(const DeformationFieldContinerType* deform, const std::vector<unsigned int>& reference, MoCoQualityType& mocoQ)
    {
        try
        {
            size_t row = deform[0].rows();
            std::vector<size_t> cols = deform[0].cols();

            long long r;

            for ( r=0; r<(long long)row; r++ )
            {
                size_t ref = reference[r];
                mocoQ[r][ref].first = 0;

                long long c;

                ValueType sumDeform = 0;
                for ( c=0; c<(long long)cols[r]; c++ )
                {
                    sumDeform += mocoQ[r][c].first;
                }

                if ( !moco_quality_determined_in_ref_picking_ || std::abs(sumDeform)<0.001 )
                {
                    {
                        #pragma omp parallel default(none) private(c) shared(cols, r, deform, mocoQ, ref) if (D==2)
                        {
                            TransformationDeformationFieldType deformTransform;
                            DeformationFieldType* deformField[D];

                            hoNDArray<CoordType> jac2D;
                            CoordType meanDeform, maxDeform, meanLogJac2D, maxLogJac2D;

                            #pragma omp for 
                            for ( c=0; c<(long long)cols[r]; c++ )
                            {
                                if ( c == ref ) continue;

                                unsigned int ii;
                                for ( ii=0; ii<D; ii++ )
                                {
                                    deformField[ii] = const_cast<DeformationFieldType*>( &(deform[ii](r, c)) );
                                }

                                deformTransform.jacobianPosition(jac2D, deformField, 1);
                                deformTransform.analyzeJacobianAndDeformation(jac2D, deformField, meanDeform, maxDeform, meanLogJac2D, maxLogJac2D, 1);

                                // mocoQ[r][c].first = meanLogJac2D;
                                mocoQ[r][c].first = meanDeform;
                            }
                        }
                    }
                }

                std::sort(mocoQ[r].begin(), mocoQ[r].end(), Gadgetron::compObj() );
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualityDeformation(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool GtMoCoAveraging<ValueType, CoordType, D>::computeMoCoQuality()
    {
        try
        {
            if ( moco_quality_strategy_ == "SSD" )
            {
                GADGET_CHECK_RETURN_FALSE(computeMoCoQualitySSD(register_.warped_container_, reference_, moco_quality_measure_));
            }
            else if ( moco_quality_strategy_ == "Deformation" )
            {
                GADGET_CHECK_RETURN_FALSE(computeMoCoQualityDeformation(register_.deformation_field_, reference_, moco_quality_measure_));
            }
            else
            {
                GERROR_STREAM("Incorrect moco quality measurement mode : " << moco_quality_strategy_);
                return false;
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GtMoCoAveraging<ValueType, CoordType, D>::computeMoCoQuality(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    void GtMoCoAveraging<ValueType, CoordType, D>::print(std::ostream& os) const
    {
        using namespace std;

        os << "-------------- Gt MoCo Averaging -------------" << endl;

        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;

        elemTypeName = std::string(typeid(CoordType).name());
        os << "Transformation coordinate data type is : " << elemTypeName << std::endl;

        os << "Whether to perform cross row registration is : " << cross_row_reg_ << std::endl;
        os << "Method to compute moco quality is : " << moco_quality_strategy_ << std::endl;
        os << "Determine moco quality in ref picking stage is : " << moco_quality_determined_in_ref_picking_ << std::endl;
        os << "Averaging percentage is : " << percentage_kept_for_averaging_ << std::endl;
        os << "Whether to perform soft averaging is : " << soft_averaging_ << std::endl;
        os << "--------------------------------------------------" << endl;
        os << "Info of register is : " << endl;
        os << "--------------------------------------------------" << endl;
        register_.print(os);
        os << "--------------------------------------------------" << endl;
    }

    // ------------------------------------------------------------
    // Instantiation
    // ------------------------------------------------------------

    template class EXPORTGTTOOLBOXCORE GtMoCoAveraging< float, double, 2 >;
    template class EXPORTGTTOOLBOXCORE GtMoCoAveraging< float, double, 3 >;
}
