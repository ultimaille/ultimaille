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

    inline int Tetrahedra::cell_type() const {
        return 0;
    }

    inline int Tetrahedra::nverts_per_cell() const {
        return 4;
    }

    inline int Tetrahedra::nfacets_per_cell() const {
        return 4;
    }

    inline int Tetrahedra::ncells() const {
        return cells.size()/4;
    }

    inline int Tetrahedra::nfacets() const {
        return ncells()*4;
    }

    inline int Tetrahedra::ncorners() const {
        return ncells()*4;
    }

    int Tetrahedra::create_tets(const int n) {
        cells.resize(cells.size()+n*4);
        resize_attrs();
        return ncells()-n;
    }

    inline int Tetrahedra::facet_size(const int c, const int lf) const {
        return 3;
    }

    int Tetrahedra::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<4 && lv>=0 && lv<3);
        const int facet_vertex[4][3] = {{1,3,2}, {0,2,3}, {3,1,0}, {0,1,2}};
        return vert(c, facet_vertex[lf][lv]);
    }

    inline int Tetrahedra::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<4);
        return c*4 + lf;
    }

    inline int Tetrahedra::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<4);
        return c*4 + lc;
    }

//  inline int Tetrahedra::vert(const int c, const int lv) const {
//      return cells[corner(c, lv)];
//  }

//  inline int &Tetrahedra::vert(const int c, const int lv) {
//      return cells[corner(c, lv)];
//  }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Hexahedra::cell_type() const {
        return 1;
    }

    inline int Hexahedra::nverts_per_cell() const {
        return 8;
    }

    inline int Hexahedra::nfacets_per_cell() const {
        return 6;
    }

    inline int Hexahedra::ncells() const {
        return cells.size()/8;
    }

    inline int Hexahedra::nfacets() const {
        return ncells()*6;
    }

    inline int Hexahedra::ncorners() const {
        return ncells()*8;
    }

    int Hexahedra::create_hexa(const int n) {
        cells.resize(cells.size()+n*8);
        resize_attrs();
        return ncells()-n;
    }

    inline int Hexahedra::facet_size(const int c, const int lf) const {
        return 4;
    }

    int Hexahedra::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<6 && lv>=0 && lv<4);
        const int facet_vertex[6][4] = {{0,2,6,4}, {3,1,5,7}, {1,0,4,5}, {2,3,7,6}, {1,3,2,0}, {4,6,7,5}};
        return vert(c, facet_vertex[lf][lv]);
    }

    inline int Hexahedra::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<6);
        return c*6 + lf;
    }

    inline int Hexahedra::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<8);
        return c*8 + lc;
    }

//  inline int Hexahedra::vert(const int c, const int lv) const {
//      return cells[corner(c, lv)];
//  }

//  inline int &Hexahedra::vert(const int c, const int lv) {
//      return cells[corner(c, lv)];
//  }


    VolumeConnectivity::VolumeConnectivity(const Volume &p_m) : m(p_m), adjacent(m.nfacets(), -1) {
        std::vector<int> c2c(m.ncorners(), -1);
        std::vector<int> v2c(m.nverts(),   -1);

        // step 1: compute vertex-to-cell corner map
        for (int c=0; c<m.ncells(); c++)
            for (int lv=0; lv<m.nverts_per_cell(); lv++)
                v2c[m.vert(c, lv)] = m.corner(c, lv);

        // step 2: chain cell corners around vertices
        for (int c=0; c<m.ncells(); c++) {
            for (int lv=0; lv<m.nverts_per_cell(); lv++) {
                int c = m.corner(c, lv);
                int v = m.vert(c, lv);
                c2c[c] = v2c[v];
                v2c[v] = c;
            }
        }

        // step 3: connect cells
        for (int c1=0; c1<m.ncells(); c1++) {
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                if (adjacent[m.facet(c1, lf1)]>=0) continue;
//              int c = m.corner(c1, lv);
//              int v = m.vert(c1, lv);
            }
        }

/*
        if(nb() == 0) {
            return;
        }
        cell_facets_.resize_store(nb() * 4);
        for(index_t f=0; f<cell_facets_.nb(); ++f) {
            cell_facets_.set_adjacent_cell(f,NO_CELL);
        }

        GEO::vector<index_t> next_tet_corner_around_vertex(
            nb() * 4, NO_CORNER
        );
        GEO::vector<index_t> v2c(vertices_.nb(), NO_CORNER);

        // Step 1: chain tet corners around vertices and compute v2c
        for(index_t t = 0; t < nb(); ++t) {
            for(index_t lv = 0; lv < 4; ++lv) {
                index_t v = vertex(t, lv);
                next_tet_corner_around_vertex[4 * t + lv] = v2c[v];
                v2c[v] = 4 * t + lv;
            }
        }

        // Step 2: connect tets
        for(index_t t1 = 0; t1 < nb(); ++t1) {
            for(index_t lf1 = 0; lf1 < 4; ++lf1) {
                if(adjacent(t1, lf1) == NO_CELL) {
                    index_t v1 = facet_vertex(t1, lf1, 0);
                    index_t v2 = facet_vertex(t1, lf1, 1);
                    index_t v3 = facet_vertex(t1, lf1, 2);
                    for(
                        index_t c2 = v2c[v1]; c2 != NO_CORNER;
                        c2 = next_tet_corner_around_vertex[c2]
                    ) {
                        index_t t2 = c2/4;
                        index_t lf2 = find_tet_facet(t2, v3, v2, v1);
                        if(lf2 != NO_FACET) {
                            set_adjacent(t1, lf1, t2);
                            set_adjacent(t2, lf2, t1);
                            break;
                        }
                    }
                }
            }
        }
*/


//      c2f.resize(nbc, -1);
//      c2c.resize(nbc, -1);
//      v2c.resize(m.nverts(), -1);
    }
}

