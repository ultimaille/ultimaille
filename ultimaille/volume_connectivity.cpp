#include <iostream>
#include <algorithm>
#include <cassert>

#include "volume.h"
#include "volume_connectivity.h"
#include "attributes.h"

namespace UM {
    void compute_corner_to_corner_map(const Volume &m, std::vector<int> &v2c, std::vector<int> &c2c);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool are_facets_adjacent(const Volume &m, int c1, int c2, int lf1, int lf2) {
        int n = m.facet_size(lf1);
        if (n!=m.facet_size(lf2)) return false;

        for (int i=0; i<n; i++) {
            bool found = true;
            for (int j=0; found && j<n; j++)
                found = (m.facet_vert(c1, lf1, (i+j)%n) == m.facet_vert(c2, lf2, n-j-1));
            if (found) return true;
        }

        return false;
    }

    OppositeFacet::OppositeFacet(const Volume &m) : m(m) {
        reset();
    }

    void OppositeFacet::reset() {
        std::cerr << "Reset!" << std::endl;
        adjacent = std::vector(m.nfacets(), -1);

        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(m, v2c, c2c);

        for (int c1=0; c1<m.ncells(); c1++)
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                int f1 = m.facet(c1, lf1);
                if (adjacent[f1]>=0) continue;
                int crn1 = v2c[m.facet_vert(c1, lf1, 0)];
                int crn2 = crn1;
                do {
                    int c2 = crn2 / m.nverts_per_cell();
                    if (c2!=c1) {
                        for (int lf2=0; lf2<m.nfacets_per_cell(); lf2++) {
                            if (!are_facets_adjacent(m, c1, c2, lf1, lf2)) continue;
                            int f2 = m.facet(c2, lf2);
                            adjacent[f1] = f2;
                            adjacent[f2] = f1;
                        }
                    }
                    crn2 = c2c[crn2];
                } while (crn2!=crn1 && adjacent[f1]<0);
            }
    }

    int & OppositeFacet::operator[](const int i)       { assert(i>=0 && i<m.nfacets()); return adjacent[i]; }
    int   OppositeFacet::operator[](const int i) const { assert(i>=0 && i<m.nfacets()); return adjacent[i]; }

    int OppositeFacet::opposite_c(const int he) const {
        assert(he>=0 && he<m.heh.nhalfedges());
        int hfacet = m.heh.facet(he);
        assert((int)adjacent.size()>hfacet);
        int opp_hfacet = adjacent[hfacet];
        if (opp_hfacet<0) return -1;
        int opp_cell = m.cell_from_facet(opp_hfacet);
        return m.heh.halfedge_from_verts(opp_cell, m.heh.to(he), m.heh.from(he));
    }

    std::vector<int> OppositeFacet::halfedges_around_edge(const int he) const {
        std::vector<int> result;
        int around_e_cir = he;
        do { // rewind if boundary
            result.push_back(around_e_cir);
            if (opposite_c(around_e_cir)<0) break;
            around_e_cir = m.heh.opposite_f(opposite_c(around_e_cir));
        } while (around_e_cir != he);

        if (opposite_c(around_e_cir)>=0 && around_e_cir == he)
            return result;

        result.clear(); // iterate forward if the edge is on border
        do {
            result.push_back(around_e_cir);
            around_e_cir = opposite_c(m.heh.opposite_f(around_e_cir));
        } while (around_e_cir>=0);
        return result;
    }

    halfedge_around_edge_iter::halfedge_around_edge_iter(const OppositeFacet &of, const int he) : of(of), heh(of.m.heh) {
        assert(he>=0 && he<heh.nhalfedges());
        int cur = he;
        do {
            int oppc = of.opposite_c(cur);
            if (oppc<0) {
                start = cur;
                return;
            }
            cur = heh.opposite_f(oppc);
        } while (cur != he);
        start = finish = he;
    }


}

