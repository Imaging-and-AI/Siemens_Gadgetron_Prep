/** \file jointT1T2RecoveryOperator.h
    \brief Implement joint T1/T2 SASHA recovery operator
    \author     Kelvin Chow
*/

#pragma once

#include <gadgetron/curveFittingOperator.h>
#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron {

    // A0*(1-exp(-d/A1)*(1 - exp(-e/A2) + exp(-s/A1)*exp(-e/A2)))
    // d/dz x*(1-e^(-d/y)*(1 - e^(-f/z) + e^(-s/y)*e^(-f/z)))

    template <class ARRAY> class jointSashaT1T2RecoveryOperator : public curveFittingOperator<ARRAY>
    {
    public:

        typedef curveFittingOperator<ARRAY> BaseClass;
        typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
        typedef typename BaseClass::REAL REAL;

        jointSashaT1T2RecoveryOperator();
        virtual ~jointSashaT1T2RecoveryOperator();

        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad);
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y);

    protected:
    };

    template <class ARRAY>
    jointSashaT1T2RecoveryOperator<ARRAY>::jointSashaT1T2RecoveryOperator() : BaseClass()
    {
    }

    template <class ARRAY>
    jointSashaT1T2RecoveryOperator<ARRAY>::~jointSashaT1T2RecoveryOperator()
    {
    }

    template <class ARRAY>
    void jointSashaT1T2RecoveryOperator<ARRAY>::gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad)
    {
        GADGET_THROW("jointSashaT1T2RecoveryOperator<ARRAY>::gradient(...) is not yet implemented ...");

        // S = A*(1-(1-(1-exp(-x/B))*exp(-y/C))*exp(-z/B)), where:
        //   x = TS,  B = T1
        //   y = TE,  C = T2
        //   z = del

        // = d/dA(A*(1-(1-(1-exp(-x/B))*exp(-y/C))*exp(-z/B)))

        // // Ignore changes in reciprocal of T1/T2 below a certain threshold
        // int sign_b1 = boost::math::sign(b[1]);
        // ELEMENT_TYPE rb1 = 1.0 / ( (std::abs(b[1])<FLT_EPSILON) ? sign_b1*FLT_EPSILON : b[1] );

        // int sign_b2 = boost::math::sign(b[2]);
        // ELEMENT_TYPE rb2 = 1.0 / ( (std::abs(b[2])<FLT_EPSILON) ? sign_b2*FLT_EPSILON : b[2] );

        // // Simplify common terms
        // ELEMENT_TYPE Ed = exp(-1 * del * rb1);
        // ELEMENT_TYPE E1 = exp(-1 * ts[ii] * rb1);
        // ELEMENT_TYPE E2 = exp(-1 * te[ii] * rb2);



        // // Apply signal model
        // size_t ii;
        // for (ii=0; ii<num; ii++)
        // {
        //     // Simplify common terms
        //     ELEMENT_TYPE E1 = exp(-1 * ts[ii] * rb1);
        //     ELEMENT_TYPE E2 = exp(-1 * te[ii] * rb2);

        //     y[ii] = b[0] * (1 - Ed*(1 - E2 + E1*E2));
        // }

        // try
        // {
        //     // b[0] is scale factor
        //     // b[1] is T1
        //     // b[2] is T2

        //     size_t num = y.size();
        //     if (x.size() != num*2+1)
        //     {
        //         GADGET_THROW("Size mismatch between x and y in jointSashaT1T2RecoveryOperator");
        //     }

        //     // Extract values from x array
        //     ARRAY ts, te;
        //     ELEMENT_TYPE del;

        //     for (size_t i=0; i<num; i++)
        //     {
        //         ts[i] = x[i];
        //         te[i] = x[i+num];
        //     }
        //     del = x[2*num];

        //     // Ignore changes in reciprocal of T1/T2 below a certain threshold
        //     int sign_b1 = boost::math::sign(b[1]);
        //     ELEMENT_TYPE rb1 = 1.0 / ( (std::abs(b[1])<FLT_EPSILON) ? sign_b1*FLT_EPSILON : b[1] );

        //     int sign_b2 = boost::math::sign(b[2]);
        //     ELEMENT_TYPE rb2 = 1.0 / ( (std::abs(b[2])<FLT_EPSILON) ? sign_b2*FLT_EPSILON : b[2] );

        //     // Simplify common terms
        //     ELEMENT_TYPE Ed = exp(-1 * del * rb1);


        //     size_t num = b.size();
        //     if(grad.size()!=num) grad.resize(num, 0);

        //     int sign_b2 = boost::math::sign(b[2]);
        //     ELEMENT_TYPE rb = 1.0 / ((std::abs(b[2])<FLT_EPSILON) ? sign_b2*FLT_EPSILON : b[2]);

        //     ELEMENT_TYPE val = std::exp(-1 * xi * rb);
        //     grad[0] = 1;
        //     grad[1] = -val;
        //     grad[2] = -1 * b[1] * val * xi * rb * rb;
        // }
        // catch(...)
        // {
        //     GADGET_THROW("Errors happened in jointSashaT1T2RecoveryOperator<ARRAY>::gradient(...) ... ");
        // }
    }

    template <class ARRAY>
    void jointSashaT1T2RecoveryOperator<ARRAY>::magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y)
    {
        // b[0] is scale factor
        // b[1] is T1
        // b[2] is T2

        // x is a vector containing TS, T2p, time-t2p-to-center, and T2p duration, of size 2N+2, where N is the size of y.
        // For the 'i'th measurement in y:
        //   x[i]     is the sat recovery time
        //   x[i+N]   is the T2p time
        //   x[end-1] is the timeT2pToCenter
        //   x[end]   is the t2pRfDuration
        ARRAY ts, te;
        ELEMENT_TYPE timeT2pToCenter, t2pRfDuration;

        size_t num = (x.size()-2)/2;

        if (y.size()!=num)
        {
            y.resize(num, 0);
        }
        ts.resize(num, 0);
        te.resize(num, 0);

        for (size_t i=0; i<num; i++)
        {
            ts[i] = x[i];
            te[i] = x[i+num];
        }
        t2pRfDuration   = x[2*num];
        timeT2pToCenter = x[2*num+1];

        for (size_t i=0; i<num; i++)
        {
            // Correct TS values for time to center
            ts[i] -= timeT2pToCenter;

            // For T2p images, subtract the duration of the T2p itself, as no T1 recovery is happening during this time (mostly)
            if (te[i] != 0)
            {
                ts[i] -= t2pRfDuration;
            }
        }

        // Ignore changes in reciprocal of T1/T2 below a certain threshold
        int sign_b1 = boost::math::sign(b[1]);
        ELEMENT_TYPE rb1 = 1.0 / ( (std::abs(b[1])<FLT_EPSILON) ? sign_b1*FLT_EPSILON : b[1] );

        int sign_b2 = boost::math::sign(b[2]);
        ELEMENT_TYPE rb2 = 1.0 / ( (std::abs(b[2])<FLT_EPSILON) ? sign_b2*FLT_EPSILON : b[2] );

        // Simplify common terms
        ELEMENT_TYPE Ed = exp(-1 * timeT2pToCenter * rb1);

        // Apply signal model
        size_t ii;
        for (ii=0; ii<num; ii++)
        {
            // Simplify common terms
            ELEMENT_TYPE E1 = exp(-1 * ts[ii] * rb1);
            ELEMENT_TYPE E2 = exp(-1 * te[ii] * rb2);

            y[ii] = b[0] * (1 - Ed*(1 - E2 + E1*E2));
        }

    }
}
