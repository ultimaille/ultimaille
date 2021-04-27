#ifndef __SURFACE_CONNECTIVITY_H__
#define __SURFACE_CONNECTIVITY_H__
#include <vector>
#include <memory>
#include "surface.h"
#include "syntactic-sugar/assert.h"

namespace UM {
    struct SurfaceConnectivity { // half-edge-like connectivity interface
        SurfaceConnectivity(const Surface &p_m);

        vec3 geom(const int halfedge) const;
        int facet(const int halfedge) const;
        int  from(const int halfedge) const;
        int    to(const int halfedge) const;
        int  prev(const int halfedge) const;
        int  next(const int halfedge) const;
        int opposite(const int halfedge) const;
        bool is_boundary_vert(const int v) const;
        int next_around_vertex(const int halfedge) const;
        int prev_around_vertex(const int halfedge) const;

        void reset();
        const Surface &m;
        std::vector<int> v2c; // vertex to corner
        std::vector<int> c2f; // corner to facet
        std::vector<int> c2c; // corner to next corner sharing the same vertex (unordered)
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // these implementations are here and not in the .cpp because all inline functions must be available in all translation units //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline vec3 SurfaceConnectivity::geom(const int halfedge) const {
        return m.points[to(halfedge)] - m.points[from(halfedge)];
    }

    inline int SurfaceConnectivity::facet(const int halfedge) const {
        return c2f[halfedge];
    }

    inline int SurfaceConnectivity::from(const int halfedge) const {
        int fi = c2f[halfedge];
        int lv = halfedge - m.corner(fi, 0);
        return m.vert(fi, lv);
    }

    inline int SurfaceConnectivity::to(const int halfedge) const {
        int fi = c2f[halfedge];
        int lv = halfedge - m.corner(fi, 0);
        int n = m.facet_size(fi);
        return m.vert(fi, (lv+1)%n);
    }

    inline int SurfaceConnectivity::next(const int halfedge) const {
        int fi = c2f[halfedge];
        int lv = halfedge - m.corner(fi, 0);
        int n = m.facet_size(fi);
        return m.corner(fi, (lv+1)%n);
    }

    inline int SurfaceConnectivity::prev(const int halfedge) const {
        int fi = c2f[halfedge];
        int lv = halfedge - m.corner(fi, 0);
        int n = m.facet_size(fi);
        return m.corner(fi, (lv-1+n)%n);
    }
}

#endif //__SURFACE_CONNECTIVITY_H__

