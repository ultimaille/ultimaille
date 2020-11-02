#include <iostream>
#include <algorithm>
#include <cassert>
#include "volume.h"
#include "attributes.h"

namespace UM {
    void Volume::resize_attrs() {
        for (auto &wp : attr_cells)   if (auto spt = wp.lock())
            spt->resize(ncells());
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->resize(nfacets());
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->resize(ncorners());
    }

    void Volume::compress_attrs(const std::vector<bool> &cells_to_kill) {
    /*
        assert(facets_to_kill.size()==(size_t)nfacets());
        std::vector<int>  facets_old2new(nfacets(),  -1);
        std::vector<int> corners_old2new(ncorners(), -1);

        int new_nb_facets  = 0;
        int new_nb_corners = 0;

        for (int f=0; f<nfacets(); f++) {
            if (facets_to_kill[f]) continue;
            for (int lv=0; lv<facet_size(f); lv++)
                corners_old2new[facet_corner(f, lv)] = new_nb_corners++;
            facets_old2new[f] = new_nb_facets++;
        }
        std::cerr << "compressing facet attributes\n";
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->compress(facets_old2new);
        std::cerr << "compressing corner attributes\n";
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->compress(corners_old2new);
       */
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Tetrahedra::create_tets(const int n) {
        cells.resize(cells.size()+n*4);
        resize_attrs();
        return ncells()-n;
    }

    int Hexahedra::create_hexa(const int n) {
        cells.resize(cells.size()+n*8);
        resize_attrs();
        return ncells()-n;
    }

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
        std::vector<int> c2c(m.ncorners(), -1);
        std::vector<int> v2c(m.nverts(),   -1);

        // step 1: compute vertex-to-cell corner map
        for (int c=0; c<m.ncells(); c++)
            for (int lv=0; lv<m.nverts_per_cell(); lv++)
                v2c[m.vert(c, lv)] = m.corner(c, lv);

        // step 2: chain cell corners around vertices
        for (int c=0; c<m.ncells(); c++)
            for (int lv=0; lv<m.nverts_per_cell(); lv++) {
                int crn = m.corner(c, lv);
                int v = m.vert(c, lv);
                c2c[crn] = v2c[v];
                v2c[v] = crn;
            }

        // step 3: connect cells
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
}

