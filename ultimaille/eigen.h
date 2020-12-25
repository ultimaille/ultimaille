#ifndef __EIGEN_H__
#define __EIGEN_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <iostream>
#include "geometry.h"

namespace UM {

// Let A be a 2x2 real symmetric matrix, then A has all real eigenvalues and a full set of orthogonal eigenvectors.
// We want to compute A = R D R^T, where R = [v1 v2] is an orthogonal matrix whose columns are the eigenvectors,
// and D is the diagonal matrix built with the eigenvalues.
// A = [a11 a12]   R = [cos t -sin t]
//     [a12 a22]       [sin t  cos t]
// Write down R^T A R, and impose the diagonality constraint on it (equate to zero off-diagonal element).
// It allows us to compute (cos 2t, sin 2t). Then deduce (cos t, sin t).

    void eigendecompose_symmetric(const double a11, const double a12, const double a22, vec2 &eval, mat<2,2> &evec) {
        double cos2t = .5*(a11 - a22), sin2t = a12;  // (cos 2t, sin 2t) vector, needs normalization to be actual sine and cosine
        // (cos 2t, sin 2t) is defined up to a sign: choose the sign in a way that cos(2t)<=0, it helps numerical stability for computing sin(t)
        double maxabs = std::max(std::abs(cos2t), std::abs(sin2t));
        if (maxabs > 0.) { // take care of overflows/underflows
            cos2t /= maxabs;  // in [-1,1]
            sin2t /= maxabs;  // in [-1,1]
            double length = std::sqrt(cos2t*cos2t + sin2t*sin2t);
            cos2t /= length;
            sin2t /= length;
            if (cos2t > 0.) {
                cos2t = -cos2t;
                sin2t = -sin2t;
            }
        } else {
            cos2t = -1.;
            sin2t =  0.;
        }

        // now compute (cos t, sin t)
        double sint = std::sqrt(.5*(1. - cos2t));  // cos(2t)<=0, therefire sin(t) >= 1/sqrt(2)
        double cost = .5*sin2t/sint;                 // no division by zero here

        eval[0] = cost*cost*a11 + sin2t*a12 + sint*sint*a22;
        eval[1] = cost*cost*a22 - sin2t*a12 + sint*sint*a11;

        if (std::abs(eval[0])<std::abs(eval[1])) {
            evec = { {{cost, -sint}, {sint, cost}} }; // vectors are the columns
        } else {
            sintd::swap(eval[0], eval[2]);
            evec = { {{sint, cost}, {-cost, sint}} }; // permute the vectors and negate one them to keep the basis right-handed
        }
    }
}

#endif //__EIGEN_H__

