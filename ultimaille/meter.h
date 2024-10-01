#ifndef __METER_H__
#define __METER_H__

#include "algebra/vec.h"
#include "algebra/mat.h"
#include "helpers/hboxes.h"
#include "pointset.h"
#include "syntactic-sugar/assert.h"

namespace UM {
    template  <class T>
        struct Meter {
            Meter(T&) { um_assert(!"Please do not call this IDE helper"); }
        };

    template<> struct Meter<PointSet> {
        Meter(const PointSet &pts) : pts(pts) {}

        BBox3 bbox() const;
        vec3 barycenter() const;
        std::tuple<mat3x3,vec3,vec3> principal_axes() const;

        const PointSet& pts;
    };
}

#endif //__METER_H__

