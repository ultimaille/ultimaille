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

    int Volume::nverts() const {
        return points.size();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Tetrahedra::cell_type() const {
        return 0;
    }

    int Tetrahedra::nverts_per_cell() const {
        return 4;
    }

    int Tetrahedra::nfacets_per_cell() const {
        return 4;
    }

    int Tetrahedra::ncells() const {
        return cells.size()/4;
    }

    int Tetrahedra::nfacets() const {
        return ncells()*4;
    }

    int Tetrahedra::ncorners() const {
        return ncells()*4;
    }

    int Tetrahedra::create_tets(const int n) {
        cells.resize(cells.size()+n*4);
        resize_attrs();
        return ncells()-n;
    }

    int Tetrahedra::vert(const int c, const int lv) const {
        assert(c>=0 && c<ncells() && lv>=0 && lv<4);
        return cells[c*4 + lv];
    }

    int &Tetrahedra::vert(const int c, const int lv) {
        assert(c>=0 && c<ncells() && lv>=0 && lv<4);
        return cells[c*4 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Hexahedra::cell_type() const {
        return 1;
    }

    int Hexahedra::nverts_per_cell() const {
        return 8;
    }

    int Hexahedra::nfacets_per_cell() const {
        return 6;
    }

    int Hexahedra::ncells() const {
        return cells.size()/8;
    }

    int Hexahedra::nfacets() const {
        return ncells()*6;
    }

    int Hexahedra::ncorners() const {
        return ncells()*8;
    }

    int Hexahedra::create_hexa(const int n) {
        cells.resize(cells.size()+n*8);
        resize_attrs();
        return ncells()-n;
    }

    int Hexahedra::vert(const int c, const int lv) const {
        assert(c>=0 && c<ncells() && lv>=0 && lv<8);
        return cells[c*8 + lv];
    }

    int &Hexahedra::vert(const int c, const int lv) {
        assert(c>=0 && c<ncells() && lv>=0 && lv<8);
        return cells[c*8 + lv];
    }
}

