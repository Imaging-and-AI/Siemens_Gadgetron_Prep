/** \file jointT1T2T1pRecoveryOperator.h
    \brief Implement joint T1/T2/T1p SASHA recovery operator
    \author     Kelvin Chow
*/

#pragma once

#include <curveFittingOperator.h>
#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron {

    template <class ARRAY> class jointSashaT1T2T1pRecoveryOperator : public curveFittingOperator<ARRAY>
    {
    public:

        typedef curveFittingOperator<ARRAY> BaseClass;
        typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
        typedef typename BaseClass::REAL REAL;

        jointSashaT1T2T1pRecoveryOperator();
        virtual ~jointSashaT1T2T1pRecoveryOperator();

        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad);
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y);

    protected:
    };

    template <class ARRAY>
    jointSashaT1T2T1pRecoveryOperator<ARRAY>::jointSashaT1T2T1pRecoveryOperator() : BaseClass()
    {
    }

    template <class ARRAY>
    jointSashaT1T2T1pRecoveryOperator<ARRAY>::~jointSashaT1T2T1pRecoveryOperator()
    {
    }

    template <class ARRAY>
    void jointSashaT1T2T1pRecoveryOperator<ARRAY>::gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad)
    {
        GADGET_THROW("jointSashaT1T2T1pRecoveryOperator<ARRAY>::gradient(...) is not yet implemented ...");
    }

    template <class ARRAY>
    void jointSashaT1T2T1pRecoveryOperator<ARRAY>::magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y)
    {
        // b[0] is scale factor
        // b[1] is T1
        // b[2] is T2
        // b[3] is T1p

        // GADGET_THROW("jointSashaT1T2T1pRecoveryOperator<ARRAY>::magnitude(...) is not yet implemented ...");

        // x is a vector containing TS, TE, TSL, time-t2p-to-center, and T2p duration, of size 4N+1, where N is the size of y.
        // For the 'i'th measurement in y:
        //   x[i]     is the sat recovery time
        //   x[i+ N]  is the T2p time
        //   x[i+2N]  is the T1p time
        //   x[i+3N]  is the t2pRfDuration
        //   x[end]   is the timeT2pToCenter

        ARRAY ts, te, tsl, t2pRfDuration;
        ELEMENT_TYPE timeT2pToCenter;

        size_t num = (x.size()-1)/4;

        if (y.size()!=num)
        {
            y.resize(num, 0);
        }
        ts.resize(           num, 0);
        te.resize(           num, 0);
        tsl.resize(          num, 0);
        t2pRfDuration.resize(num, 0);

        for (size_t i=0; i<num; i++)
        {
            ts[           i] = x[i];
            te[           i] = x[i+num];
            tsl[          i] = x[i+num*2];
            t2pRfDuration[i] = x[i+num*3];
        }
        timeT2pToCenter = x[4*num+1];

        for (size_t i=0; i<num; i++)
        {
            // Correct TS values for time to center
            ts[i] -= timeT2pToCenter;

            // Subtract the duration of the T2p itself, as no T1 recovery is happening during this time (mostly)
            // (This value is 0 for images with T2p or T1p)
            ts[i] -= t2pRfDuration[i];
        }

        // Ignore changes in reciprocal of T1/T2 below a certain threshold
        int sign_b1 = boost::math::sign(b[1]);
        ELEMENT_TYPE rb1 = 1.0 / ( (std::abs(b[1])<FLT_EPSILON) ? sign_b1*FLT_EPSILON : b[1] );

        int sign_b2 = boost::math::sign(b[2]);
        ELEMENT_TYPE rb2 = 1.0 / ( (std::abs(b[2])<FLT_EPSILON) ? sign_b2*FLT_EPSILON : b[2] );

        int sign_b1p = boost::math::sign(b[3]);
        ELEMENT_TYPE rb1p = 1.0 / ( (std::abs(b[3])<FLT_EPSILON) ? sign_b1p*FLT_EPSILON : b[3] );

        // Simplify common terms
        ELEMENT_TYPE Ed = exp(-1 * timeT2pToCenter * rb1);

        // Apply signal model
        size_t ii;
        for (ii=0; ii<num; ii++)
        {
            // Simplify common terms
            ELEMENT_TYPE E1  = exp(-1 * ts[ ii] * rb1);
            ELEMENT_TYPE E2  = exp(-1 * te[ ii] * rb2);
            ELEMENT_TYPE E1p = exp(-1 * tsl[ii] * rb1p);

            y[ii] = b[0] * (1 - Ed*(1 - E2*E1p + E1*E2*E1p));
        }

    }
}
