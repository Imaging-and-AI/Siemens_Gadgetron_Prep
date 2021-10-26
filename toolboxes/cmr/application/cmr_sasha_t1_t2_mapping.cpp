/** \file   cmr_t1_mapping.cpp
    \brief  Implement CMR T1 mapping for 2D acquisition
    \author Hui Xue
*/


#include "cmr_sasha_t1_t2_mapping.h"
#include <log.h>

#include <hoNDArray_reductions.h>
#include <hoNDArray_elemwise.h>
#include <hoNDArray_math.h>

#include <simplexLagariaSolver.h>
#include "jointSashaT1T2RecoveryOperator.h"
#include "jointSashaT1T2T1pRecoveryOperator.h"
#include <curveFittingCostFunction.h>

#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron { 

template <typename T> 
CmrSashaT1T2Mapping<T>::CmrSashaT1T2Mapping() : BaseClass()
{
    max_iter_     =  150;
    max_fun_eval_ = 1000;
    thres_fun_    = 1e-4;

    // maximal allowed values
    max_t1_value_ = 5000;
    max_t2_value_ =  500;
}

template <typename T> 
CmrSashaT1T2Mapping<T>::~CmrSashaT1T2Mapping()
{
}

template <typename T>
void CmrSashaT1T2Mapping<T>::get_initial_guess(const VectorType& ti, const VectorType& yi, VectorType& guess)
{
    if (guess.size() != this->get_num_of_paras())
    {
        guess.resize(this->get_num_of_paras(), 0);
    }

    guess[0] =  500; // A
    guess[1] = 1200; // T1
    guess[2] =   50; // T2

    // A
    if(!yi.empty()) guess[0] = *std::max_element(yi.begin(), yi.end());
}

template <typename T>
void CmrSashaT1T2Mapping<T>::compute_map(const VectorType& ti, const VectorType& yi, const VectorType& guess, VectorType& bi, T& map_v)
{
    try
    {
        bi = guess;
        map_v = 0;

        typedef Gadgetron::jointSashaT1T2RecoveryOperator<std::vector<T> > SignalType;
        typedef Gadgetron::leastSquareErrorCostFunction<  std::vector<T> > CostType;

        // define solver
        Gadgetron::simplexLagariaSolver< VectorType, SignalType, CostType > solver;

        // define signal model
        SignalType multi_sasha;

        // define cost function
        CostType lse;

        solver.signal_model_ = &multi_sasha;
        solver.cf_           = &lse;

        solver.max_iter_     = max_iter_;
        solver.max_fun_eval_ = max_fun_eval_;
        solver.thres_fun_    = thres_fun_;

        solver.x_ = ti;
        solver.y_ = yi;

        solver.solve(bi, guess);

        // Validate calculated values
        if (bi[0] > 0 && bi[1] > 0 && bi[2] > 0)
        {
            // !!! Need to return multiple values !!!
            // T1
            if (bi[1] >= max_t1_value_)  bi[1] = hole_marking_value_;
            if (bi[1] <= min_map_value_) bi[1] = hole_marking_value_;

            if (bi[2] >= max_t2_value_)  bi[2] = hole_marking_value_;
            if (bi[2] <= min_map_value_) bi[2] = hole_marking_value_;

            map_v = bi[1];
        }
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in CmrSashaT1T2Mapping<T>::compute_map(...) ... ");
    }
}

template <typename T>
void CmrSashaT1T2Mapping<T>::compute_sd(const VectorType& ti, const VectorType& yi, const VectorType& bi, VectorType& sd, T& map_sd)
{
    GADGET_THROW("CmrSashaT1T2Mapping<T>::compute_sd(...) is not yet implemented... ");
}

template <typename T>
size_t CmrSashaT1T2Mapping<T>::get_num_of_paras() const
{
    return 3; // A, T1, T2
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class CmrSashaT1T2Mapping< float >;

// --------------------------------------------------------------

template <typename T>
CmrSashaT1T2T1pMapping<T>::CmrSashaT1T2T1pMapping() : BaseClass()
{
    max_iter_ = 150;
    max_fun_eval_ = 1000;
    thres_fun_ = 1e-4;

    // maximal allowed values
    max_t1_value_ = 5000;
    max_t2_value_ = 500;
    max_t1p_value_ = 500;
}

template <typename T>
CmrSashaT1T2T1pMapping<T>::~CmrSashaT1T2T1pMapping()
{
}

template <typename T>
void CmrSashaT1T2T1pMapping<T>::get_initial_guess(const VectorType& ti, const VectorType& yi, VectorType& guess)
{
    if (guess.size() != this->get_num_of_paras())
    {
        guess.resize(this->get_num_of_paras(), 0);
    }

    guess[0] = 500; // A
    guess[1] = 1200; // T1
    guess[2] = 50; // T2
    guess[3] = 50; // T1p

    // A
    if (!yi.empty()) guess[0] = *std::max_element(yi.begin(), yi.end());
}

template <typename T>
void CmrSashaT1T2T1pMapping<T>::compute_map(const VectorType& ti, const VectorType& yi, const VectorType& guess, VectorType& bi, T& map_v)
{
    try
    {
        bi = guess;
        map_v = 0;

        typedef Gadgetron::jointSashaT1T2T1pRecoveryOperator<std::vector<T> > SignalType;
        typedef Gadgetron::leastSquareErrorCostFunction<  std::vector<T> > CostType;

        // define solver
        Gadgetron::simplexLagariaSolver< VectorType, SignalType, CostType > solver;

        // define signal model
        SignalType multi_sasha;

        // define cost function
        CostType lse;

        solver.signal_model_ = &multi_sasha;
        solver.cf_ = &lse;

        solver.max_iter_ = max_iter_;
        solver.max_fun_eval_ = max_fun_eval_;
        solver.thres_fun_ = thres_fun_;

        solver.x_ = ti;
        solver.y_ = yi;

        solver.solve(bi, guess);

        // Validate calculated values
        if (bi[0] > 0 && bi[1] > 0 && bi[2] > 0 && bi[3] > 0)
        {
            // !!! Need to return multiple values !!!
            if (bi[1] >= max_t1_value_)  bi[1] = hole_marking_value_;
            if (bi[1] <= min_map_value_) bi[1] = hole_marking_value_;

            if (bi[2] >= max_t2_value_)  bi[2] = hole_marking_value_;
            if (bi[2] <= min_map_value_) bi[2] = hole_marking_value_;

            if (bi[3] >= max_t1p_value_) bi[3] = hole_marking_value_;
            if (bi[3] <= min_map_value_) bi[3] = hole_marking_value_;

            map_v = bi[1];
        }
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in CmrSashaT1T2T1pMapping<T>::compute_map(...) ... ");
    }
}

template <typename T>
void CmrSashaT1T2T1pMapping<T>::compute_sd(const VectorType& ti, const VectorType& yi, const VectorType& bi, VectorType& sd, T& map_sd)
{
    GADGET_THROW("CmrSashaT1T2T1pMapping<T>::compute_sd(...) is not yet implemented... ");
}

template <typename T>
size_t CmrSashaT1T2T1pMapping<T>::get_num_of_paras() const
{
    return 4; // A, T1, T2, T1p
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class CmrSashaT1T2T1pMapping< float >;

}
