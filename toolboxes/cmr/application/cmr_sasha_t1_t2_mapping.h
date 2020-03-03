/** \file   cmr_sasha_t1_t2_mapping.h
    \brief  Implement cardiac MR joint t1/T2 mapping with SASHA for 2D applications
            The input has dimension [RO E1 N S SLC]
            Temporal dimension is N
    \author Hui Xue
*/

#pragma once

#include "gadgetron_siemens_toolbox_cmr_export.h"
#include "cmr_parametric_mapping.h"
#include "cmr_multi_parametric_mapping.h"
#include "cmr_sasha_t1_t2_mapping.h"

namespace Gadgetron { 

// ======================================================================================
// SASHA Joint T1/T2
// y = b[0] * (1 - Ed*(1 - E2 + E1.*E2) );
//   where:
//     E1 = exp(-aTS / b[1]);
//     E2 = exp(-aTE / b[2]);
//     Ed = exp(-del / b[1]);

template <typename T>
class EXPORTGTTOOLBOXCMR CmrSashaT1T2Mapping : public CmrMultiParametricMapping<T>
{
public:

    typedef CmrMultiParametricMapping<T> BaseClass;
    typedef CmrSashaT1T2Mapping<T> Self;

    typedef typename BaseClass::ArrayType ArrayType;
    typedef typename BaseClass::ImageType ImageType;
    typedef typename BaseClass::ImageContinerType ImageContinerType;
    typedef typename BaseClass::VectorType VectorType;

    CmrSashaT1T2Mapping();
    virtual ~CmrSashaT1T2Mapping();

    // ======================================================================================
    /// parameter for SASHA T1/T2 mapping
    // ======================================================================================
    T max_t1_value_;
    T max_t2_value_;

    // ======================================================================================
    // perform every steps
    // ======================================================================================

    /// provide initial guess for the mapping
    virtual void get_initial_guess(const VectorType& ti, const VectorType& yi, VectorType& guess);

    /// compute map values for every parameters in bi
    virtual void compute_map(const VectorType& ti, const VectorType& yi, const VectorType& guess, VectorType& bi, T& map_v);

    /// compute SD values for every parameters in bi
    virtual void compute_sd(const VectorType& ti, const VectorType& yi, const VectorType& bi, VectorType& sd, T& map_sd);

    /// two parameters, A, T1
    virtual size_t get_num_of_paras() const;

    // ======================================================================================
    /// parameter from BaseClass
    // ======================================================================================

    using BaseClass::fill_holes_in_maps_;
    using BaseClass::max_size_of_holes_;
    using BaseClass::hole_marking_value_;
    using BaseClass::compute_SD_maps_;
    using BaseClass::mask_for_mapping_;
    using BaseClass::ti_;
    using BaseClass::data_;
    using BaseClass::map_;
    using BaseClass::para_;
    using BaseClass::sd_map_;
    using BaseClass::sd_para_;
    using BaseClass::max_iter_;
    using BaseClass::max_fun_eval_;
    using BaseClass::thres_fun_;
    using BaseClass::max_map_value_;
    using BaseClass::min_map_value_;

    using BaseClass::verbose_;
    using BaseClass::debug_folder_;
    using BaseClass::perform_timing_;
    using BaseClass::gt_timer_local_;
    using BaseClass::gt_timer_;
    using BaseClass::gt_exporter_;
};

}
