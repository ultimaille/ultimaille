#ifndef __COVARIANCE_H__
#define __COVARIANCE_H__

#include <vector>
#include "vec.h"
#include "mat.h"

namespace UM {
    struct PointSetCovariance {
        PointSetCovariance() = default;
        PointSetCovariance(const std::vector<vec3>& pts) {
            n = static_cast<int>(pts.size());
            for (const vec3 &p : pts) center += p;
            center /= static_cast<double>(n);
            for (const vec3 &p : pts)
                for (int i=0; i<3; i++)
                    for (int j=0; j<3; j++)
                        cov[i][j] += (p-center)[i]*(p-center)[j];
            cov /= static_cast<double>(n);
        }

        mat3x3 cov = {};
        vec3 center = {};
        int n = 0;
    };

    // Pebay, Philippe Pierre. Formulas for robust, one-pass parallel computation of covariances and arbitrary-order statistical moments.  doi:10.2172/1028931.
    PointSetCovariance operator+(PointSetCovariance& A, PointSetCovariance& B);
}

#endif //__COVARIANCE_H__

