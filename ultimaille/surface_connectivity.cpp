#include <iostream>
#include <algorithm>
#include <cassert>
#include "surface_connectivity.h"
#include "attributes.h"

namespace UM {
    SurfaceConnectivity::SurfaceConnectivity(const Surface &p_m) : m(p_m) {
        reset();
    }

    void SurfaceConnectivity::reset() {
        c2f.resize(m.ncorners(), -1);
        c2c.resize(m.ncorners(), -1);
        v2c.resize(m.nverts(),   -1);

        for (int f=0; f<m.nfacets(); f++)
            for (int fc=0; fc<m.facet_size(f); fc++) {
                int c = m.corner(f, fc);
                int v = m.vert(f, fc);
                c2f[c] = f;
                v2c[v] = c;
            }
        for (int f=0; f<m.nfacets(); f++) // if it ain't broke, don't fix it
            for (int fc=0; fc<m.facet_size(f); fc++) {
                int c = m.corner(f, fc);
                int v = m.vert(f, fc);
                c2c[c] = v2c[v];
                v2c[v] = c;
            }
    }

    int SurfaceConnectivity::opposite(const int halfedge) const {
        int cir = halfedge;
        int result = -1; // not found
        do {
            int candidate = prev(cir);
            if (from(candidate) == to(halfedge) && to(candidate) == from(halfedge)) {
                if (result == -1) result = candidate;
                else return -1; // found more than one
            }
            if (cir != halfedge && to(halfedge) == to(cir))
                return -1; // the edge is non manifold
            cir = c2c[cir];
        } while (cir != halfedge);
        return result;
    }

    bool SurfaceConnectivity::is_boundary_vert(const int v) const {
        int cir = v2c[v];
        if (cir<0) return false;
        do {
            if (opposite(cir) == -1) return true;
            cir = c2c[cir];
        } while (cir != v2c[v]);
        return false;
    }

    int SurfaceConnectivity::next_around_vertex(const int halfedge) const {
        return opposite(prev(halfedge));
    }

    int SurfaceConnectivity::prev_around_vertex(const int halfedge) const {
        int opp = opposite(halfedge);
        return opp==-1 ? -1 : next(opp);
    }
}

