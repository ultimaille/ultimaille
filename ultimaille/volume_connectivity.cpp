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
        adjacent = std::vector(m.nfacets(), -1); // -1 not found already, -2 more than one candidate

        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(m, v2c, c2c);

        for (int c1=0; c1<m.ncells(); c1++)
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                int f1 = m.facet(c1, lf1);
                int crn1 = v2c[m.facet_vert(c1, lf1, 0)];
                int crn2 = crn1;
                do {
                    int c2 = crn2 / m.nverts_per_cell();
                    if (c2!=c1) {
                        for (int lf2=0; lf2<m.nfacets_per_cell(); lf2++) {
                            if (!are_facets_adjacent(m, c1, c2, lf1, lf2)) continue;
                            int f2 = m.facet(c2, lf2);
                            if (adjacent[f1]==-1 && adjacent[f2]==-1 ){ // perfect match
                                adjacent[f1] = f2;
                                adjacent[f2] = f1;
                                continue;
                            }
                            // destroy previous links
                            if (adjacent[f1]>=0)  adjacent[adjacent[f1]] = -2;
                            if (adjacent[f2]>=0)  adjacent[adjacent[f2]] = -2; 

                            adjacent[f1] = -2;
                            adjacent[f2] = -2;
                        }
                    }
                    crn2 = c2c[crn2];
                } while (crn2!=crn1 && adjacent[f1]<0);
            }
        for (int c1=0; c1<m.ncells(); c1++)
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                int f1 = m.facet(c1, lf1);
                if (adjacent[f1] == -2) adjacent[f1] = -1;
            }

    }

    int OppositeFacet::operator[](const int i) const {
        assert(i>=0 && i<m.nfacets());
        return adjacent[i];
    }

}

