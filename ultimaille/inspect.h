#ifndef __INSPECT_H__
#define __INSPECT_H__

#include "algebra/vec.h"
#include "algebra/mat.h"
#include "helpers/hboxes.h"
#include "pointset.h"
#include "surface.h"
#include "syntactic-sugar/assert.h"

namespace UM {
    template  <class T>
        struct Inspect {
            Inspect(T&) { um_assert(!"Please do not call this IDE helper"); }
        };

    template<> struct Inspect<PointSet> {
        Inspect(const PointSet &pts) : pts(pts) {}

        BBox3 bbox() const;
        vec3 barycenter() const;
        std::tuple<mat3x3,vec3,vec3> principal_axes() const;

        const PointSet& pts;
    };

    template<> struct Inspect<Surface::Vertex> {
        Inspect(const Surface::Vertex v) : v(v) {}
        int valence() {
            int ret = 0;
            for ([[maybe_unused]] Surface::Halfedge cir : v.iter_halfedges()) {
                assert(cir.active());
                ret++;
            }
            return ret;
        }
        Surface::Vertex v;
    };

    template<> struct Inspect<Surface::Halfedge> {
        Inspect(const Surface::Halfedge h) : h(h) {}

        double corner_angle() const {
            const vec3 &a = h.from();
            const vec3 &b = h.to();
            const vec3 &c = h.prev().from();
            return geo::angle(b-a, c-a);
        }

        double dihedral_angle() const {
            um_assert(h.opposite().active());
            return geo::dihedral_angle(Triangle3(h.facet()).normal(), Triangle3(h.opposite().facet()).normal());
        }

        Surface::Halfedge next_on_border() const {
            for (auto& it : h.to().iter_halfedges()) if (!it.opposite().active()) return it;
            um_assert(false);
            return h;
        }

        const Surface::Halfedge h;
    };

    template<> struct Inspect<Surface> {
        Inspect(Surface &m) : m(m) {}
        bool is_manifold(bool verbose=true) const; // in the engineering sense
        bool is_disk(bool verbose=true) const;
        int nb_boundaries() const;
        int nb_connected_components() const;
        int euler_characteristic() const;
        double avg_edge_len() const;
        Surface &m;
    };

}

#endif //__INSPECT_H__

