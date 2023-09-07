#ifndef __SVD_H__
#define __SVD_H__

#include <tuple>
#include "mat.h"

namespace UM {
    // Input:  arbitrary 2x2 matrix M
    // Output: mat2x2 U,D,V s.t. M = U * D * V.transpose()
    // D is diagonal, sorted by the magnitude of singular values. The singular values are guaranteed to be positive.

    std::tuple<mat2x2,mat2x2,mat2x2> svd2x2(const mat2x2 &M);

    std::tuple<mat3x3,mat3x3,mat3x3> svd3x3(const mat3x3 &M); // TODO make a better implementation, this is a placeholder
}

#endif //__SVD_H__

