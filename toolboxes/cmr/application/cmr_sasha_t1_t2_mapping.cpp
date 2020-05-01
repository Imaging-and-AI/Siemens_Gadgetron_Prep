/** \file   cmr_t1_mapping.cpp
    \brief  Implement CMR T1 mapping for 2D acquisition
    \author Hui Xue
*/


#include "cmr_sasha_t1_t2_mapping.h"
#include <gadgetron/log.h>

#include <gadgetron/hoNDArray_reductions.h>
#include <gadgetron/hoNDArray_elemwise.h>
#include <gadgetron/hoNDArray_math.h>

#include <gadgetron/simplexLagariaSolver.h>
#include "jointSashaT1T2RecoveryOperator.h"
#include <gadgetron/curveFittingCostFunction.h>

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

    // Better guesses:
    // T1 is 
    // sortData = sort(data( (aTS ~= 1e5) & (aTE == 0) ));
    // indTest  = find(data == sortData(round(numel(sortData)/2)),1);
    // startT1 = -aTS(indTest)/log(1-data(indTest)/max(data));
    // if isnan(startT1) || ~isreal(startT1)
    // 	startT1 = 1000;
    // end

    // VectorType res_sorted(res);
    // std::sort(res_sorted.begin(), res_sorted.end());

    // T1
    // if (!ti.empty()) guess[1] = ti[ti.size() / 2];
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

        // for (size_t n = 0; n < solver.x_.size(); n++)
        // {
        //     GDEBUG_STREAM("solver.x_[" << n << "] = " << solver.x_[n]);
        // }

        // for (size_t n = 0; n < solver.y_.size(); n++)
        // {
        //     GDEBUG_STREAM("solver.y_[" << n << "] = " << solver.y_[n]);
        // }

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
    // try
    // {
    //     sd.clear();
    //     sd.resize(bi.size(), 0);

    //     map_sd = 0;

    //     typedef Gadgetron::twoParaExpRecoveryOperator< std::vector<T> > SignalType;
    //     SignalType multi_sasha;

    //     // compute fitting values
    //     VectorType y;
    //     multi_sasha.magnitude(ti, bi, y);

    //     // compute residual
    //     VectorType res(y), abs_res(y);

    //     size_t num = ti.size();
    //     size_t N = this->get_num_of_paras();

    //     size_t n;
    //     for (n = 0; n < num; n++)
    //     {
    //         res[n] = y[n] - yi[n];
    //         abs_res[n] = std::abs(res[n]);
    //     }

    //     hoNDArray<T> grad;
    //     grad.create(N, num);
    //     Gadgetron::clear(grad);

    //     VectorType gradVec(N);
    //     for (n = 0; n < num; n++)
    //     {
    //         multi_sasha.gradient(ti[n], bi, gradVec);
    //         memcpy(grad.begin() + n*N, &gradVec[0], sizeof(T)*N);
    //     }

    //     GADGET_CATCH_THROW(this->compute_sd_impl(ti, yi, bi, abs_res, grad, sd));

    //     map_sd = sd[1];
    //     if (map_sd > max_map_value_) map_sd = this->hole_marking_value_;
    // }
    // catch (...)
    // {
    //     GADGET_THROW("Exceptions happened in CmrSashaT1T2Mapping<T>::compute_map(...) ... ");
    // }
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

}
