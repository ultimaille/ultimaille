#ifndef __METER_H__
#define __METER_H__

#include "algebra/vec.h"
#include "algebra/mat.h"
#include "helpers/hboxes.h"
#include "pointset.h"
#include "surface.h"
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

    template<> struct Meter<Surface::Vertex> {
        Meter(const Surface::Vertex v) : v(v) {}
        int valence() {
            int ret = 0;
            for (Surface::Halfedge cir : v.iter_halfedges()) {
                assert(cir.active());
                ret++;
            }
            return ret;
        }
        Surface::Vertex v;
    };

    template<> struct Meter<Surface::Halfedge> {
        Meter(const Surface::Halfedge h) : h(h) {}

        double corner_angle() const {
            return geo::angle(h.to().pos() - h.from().pos(), h.prev().from().pos() - h.from().pos());
        }

        double dihedral_angle() const {
            um_assert(h.opposite().active());
            return std::acos(std::clamp(Triangle3(h.facet()).normal() * Triangle3(h.opposite().facet()).normal(), -1., 1.));
        }

        Surface::Halfedge next_on_border() const {
            for (auto& it : h.to().iter_halfedges()) if (!it.opposite().active()) return it;
            um_assert(false);
            return h;
        }

        const Surface::Halfedge h;
    };

}

#endif //__METER_H__

