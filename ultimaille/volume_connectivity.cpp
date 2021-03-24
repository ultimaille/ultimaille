#include <iostream>
#include <algorithm>
#include <cassert>

#include "volume_connectivity.h"
#include "attributes.h"

namespace UM {
    void compute_corner_to_corner_map(const Volume &m, std::vector<int> &v2c, std::vector<int> &c2c);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool are_facets_adjacent(const Volume &m, int c1, int c2, int lf1, int lf2) {
        int n = m.facet_size(c1, lf1);
        if (n!=m.facet_size(c2, lf2)) return false;

        for (int i=0; i<n; i++) {
            bool found = true;
            for (int j=0; found && j<n; j++)
                found = (m.facet_vert(c1, lf1, (i+j)%n) == m.facet_vert(c2, lf2, n-j-1));
            if (found) return true;
        }

        return false;
    }

    VolumeConnectivity::VolumeConnectivity(const Volume &p_m) : m(p_m), adjacent(m.nfacets(), -1) {
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

    // c, org and dst are global indices
    int VolumeConnectivity::halfedge_from_verts(const int c, const int org, const int dst) const {
        assert(c>=0); assert(org>=0); assert(dest>=0);

        for (int cf=0; cf<m.nfacets_per_cell(); cf++) {
            int nbv = m.facet_size(c, cf);
            for (int cfh=0; cfh<nbv; cfh++) {
                if (
                        org == m.facet_vert(c, cf, cfh)
                        &&
                        dst == m.facet_vert(c, cf, (cfh+1) % nbv)
                   )
                    return halfedge(c, cf, cfh);
            }
        }
        return -1;
    }

    int VolumeConnectivity::opposite_f(const int he) const {
        assert(he>=0);
        return halfedge_from_verts(cell(he), to(he), from(he));
    }

    int VolumeConnectivity::opposite_c(const int he) const {
        assert(he>=0);
        int hfacet = facet(he);
        assert(adjacent.size()>hfacet);
        int opp_hfacet = adjacent[hfacet];
        if (opp_hfacet<0) return -1;
        int opp_cell = m.cell_from_facet(opp_hfacet);
        return halfedge_from_verts(opp_cell, to(he), from(he));
    }
}

