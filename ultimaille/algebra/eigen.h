#ifndef __EIGEN_H__
#define __EIGEN_H__

#include <tuple>
#include "vec.h"
#include "mat.h"
#include "quaternion.h"

namespace UM {
    // Input:  symmetric matrix A
    // Output: eval - vector of eigenvalues, sorted by absolute value (decreasing), i.e. eval[0] is the largest magnitude eigenvalue
    //         evec - matrix of corresponding eigenvectors (columns) such that A = evec x diag(eval) evec^T
    //                The eigenvectors are guaranteed to be orthonormal and form a right-hand basis.

    // Analytic 2x2 eigensolver
    std::tuple<vec2,mat2x2> eigendecompose_symmetric(const mat2x2 &A);
    // Iterative 3x3 eigensolver
    std::tuple<vec3,mat3x3> eigendecompose_symmetric(const mat3x3 &A);
}

#endif //__EIGEN_H__

