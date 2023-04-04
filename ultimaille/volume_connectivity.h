#ifndef __VOLUME_CONNECTIVITY_H__
#define __VOLUME_CONNECTIVITY_H__

#include <vector>
#include <memory>
#include "syntactic-sugar/assert.h"
#include "volume.h"

namespace UM {
    struct Volume;

    struct OppositeFacet { // adjacency between half-facets
        OppositeFacet(const Volume &m);

        int & operator[](const int i);
        int   operator[](const int i) const;

        int opposite_c(const int he) const;
        [[deprecated]] std::vector<int> halfedges_around_edge(const int he) const;

        const Volume &m;
        std::vector<int> adjacent;
    };

    struct HalfEdgeHelper { // half-edge-like connectivity interface
        HalfEdgeHelper(const Volume &mesh) : m(mesh) {}

        int halfedge(const int cell, const int cell_facet, const int facet_he) const;
        int halfedge_from_verts(const int c, const int org, const int dst) const;

        int nhalfedges() const;
        int nhalfedges_per_cell() const;

        vec3          geom(const int he) const;
        int           cell(const int he) const;
        int          facet(const int he) const;
        int         corner(const int he) const;
        int     cell_facet(const int he) const;
        int  cell_halfedge(const int he) const;
        int facet_halfedge(const int he) const;
        int           from(const int he) const;
        int             to(const int he) const;
        int           prev(const int he) const;
        int           next(const int he) const;

        int opposite_f(const int he) const;
        int opposite_c(const OppositeFacet &adj, const int he) const;
        const Volume &m;
    };

    struct halfedge_around_edge_iter {
        const OppositeFacet  &of;
        const HalfEdgeHelper &heh;
        int start  = -1;
        int finish = -1;

        halfedge_around_edge_iter(const OppositeFacet &of, const int he);

        struct iterator {
            const OppositeFacet  &of;
            const HalfEdgeHelper &heh;
            int value;
            bool circ;

            void operator++() {
                value = of.opposite_c(heh.opposite_f(value));
                circ = false;
            }

            int operator*() const {
                return value;
            }

            bool operator!=(const iterator& rhs) const {
                return circ || value != rhs.value;
            }
        };

        iterator begin() const { return {of, heh, start, start==finish}; }
        iterator end()   const { return {of, heh, finish, false}; }
    };
}

#endif //__VOLUME_CONNECTIVITY_H__

