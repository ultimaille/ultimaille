#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <iostream>

#include "meter.h"
#include "algebra/covariance.h"
#include "algebra/eigen.h"

namespace UM {
    BBox3 Meter<PointSet>::bbox() const {
        BBox3 bbox;
        for (vec3 const &p : pts)
            bbox.add(p);
        return bbox;
    }

    vec3 Meter<PointSet>::barycenter() const {
        vec3 ave = {0, 0, 0};
        for (const vec3 &p : pts)
            ave += p;
        return ave / static_cast<double>(pts.size());
    }

    std::tuple<mat3x3,vec3,vec3> Meter<PointSet>::principal_axes() const {
        PointSetCovariance cov(*pts.data);
        auto [eval, evec] = eigendecompose_symmetric(cov.cov);
        if (pts.size()<4) return {mat3x3::identity(), {1.,1.,1.}, cov.center}; // If the system is under-determined, return the trivial basis
        return { evec, eval, cov.center };
    }
}

