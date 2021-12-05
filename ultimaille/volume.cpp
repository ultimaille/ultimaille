#include <iostream>
#include <algorithm>
#include <cassert>

#include "volume.h"
#include "attributes.h"

namespace UM {
    // signed volume for a tet with vertices (A,B,C,D) such that (AB, AC, AD) form a right hand basis
    inline double tet_volume(const vec3 &A, const vec3 &B, const vec3 &C, const vec3 &D) {
        return ((B-A)*cross(C-A,D-A))/6.;
    }

    double Volume::Util::cell_volume(const int c) const {
        if (m.cell_type==Volume::TETRAHEDRON)
            return tet_volume(
                    m.points[m.vert(c, 0)],
                    m.points[m.vert(c, 1)],
                    m.points[m.vert(c, 2)],
                    m.points[m.vert(c, 3)]
                    );

        const int nbf = m.nfacets_per_cell();
        const vec3 bary = m.util.bary_verts(c);
        double vol = 0;
        for (int lf=0; lf<nbf; lf++) {
            const int nbv = m.facet_size(lf);
            if (3==nbv) {
                vol += tet_volume(
                        bary,
                        m.points[m.facet_vert(c, lf, 0)],
                        m.points[m.facet_vert(c, lf, 1)],
                        m.points[m.facet_vert(c, lf, 2)]
                        );
            } else if (4==nbv) {
                for (int lv=0; lv<4; lv++) {
                    vol += tet_volume(
                            bary,
                            m.points[m.facet_vert(c, lf,  lv     )],
                            m.points[m.facet_vert(c, lf, (lv+1)%4)],
                            m.points[m.facet_vert(c, lf, (lv+2)%4)]
                            )*.5;
                }
            } else {
                um_assert(false);
            }
        }
        return vol;
    }

    vec3 Volume::Util::bary_verts(const int c) const {
        vec3 ave = {0, 0, 0};
        const int nbv = m.nverts_per_cell();
        for (int lv=0; lv<nbv; lv++)
            ave += m.points[m.vert(c, lv)];
        return ave / static_cast<double>(nbv);
    }

    vec3 Volume::Util::bary_facet(const int c, const int lf) const {
        vec3 ave = {0, 0, 0};
        const int nbv = m.facet_size(lf);
        for (int lv=0; lv<nbv; lv++)
            ave += m.points[m.facet_vert(c, lf, lv)];
        return ave / static_cast<double>(nbv);
    }

    // unit vector: weighted sum of normal of a triangle fan around the barycenter
    vec3 Volume::Util::facet_normal(const int c, const int lf) const {
        if (3==m.facet_size(lf)) {
            const vec3 &A = m.points[m.facet_vert(c, lf, 0)];
            const vec3 &B = m.points[m.facet_vert(c, lf, 1)];
            const vec3 &C = m.points[m.facet_vert(c, lf, 2)];
            return cross(B-A, C-A).normalized();
        }

        vec3 res = {0, 0, 0};
        vec3 bary = m.util.bary_facet(c, lf);
        const int nbv = m.facet_size(lf);
        for (int lv=0; lv<nbv; lv++)
            res += cross(
                    m.points[m.facet_vert(c, lf,  lv       )]-bary,
                    m.points[m.facet_vert(c, lf, (lv+1)%nbv)]-bary
                    );
        return res.normalized();
    }

    /////////////////////////////////////////////////////////////////

    void compute_corner_to_corner_map(const Volume &m, std::vector<int> &v2c, std::vector<int> &c2c);

    void Volume::resize_attrs() {
        for (auto &wp : attr_cells)   if (auto spt = wp.lock())
            spt->resize(ncells());
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->resize(nfacets());
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->resize(ncorners());
    }

    void Volume::compress_attrs(const std::vector<bool> &cells_to_kill) {
        assert(cells_to_kill.size()==(size_t)ncells());
        std::vector<int>   cells_old2new(ncells(),   -1);
        std::vector<int>  facets_old2new(nfacets(),  -1);
        std::vector<int> corners_old2new(ncorners(), -1);

        int new_nb_cells   = 0;
        int new_nb_facets  = 0;
        int new_nb_corners = 0;

        for (int c=0; c<ncells(); c++) {
            if (cells_to_kill[c]) continue;
            for (int lf=0; lf<nfacets_per_cell(); lf++)
                facets_old2new[facet(c, lf)] = new_nb_facets++;
            for (int lv=0; lv<nverts_per_cell(); lv++)
                corners_old2new[corner(c, lv)] = new_nb_corners++;
            cells_old2new[c] = new_nb_cells++;
        }

//      std::cerr << "compressing cell attributes\n";
        for (auto &wp : attr_cells)   if (auto spt = wp.lock())
            spt->compress(cells_old2new);
//      std::cerr << "compressing facet attributes\n";
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->compress(facets_old2new);
//      std::cerr << "compressing corner attributes\n";
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->compress(corners_old2new);
    }

    void Volume::delete_cells(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)ncells());
        compress_attrs(to_kill);

        int new_nb_cells   = 0;
        int new_nb_corners = 0;
        for (int c=0; c<ncells(); c++) {
            if (to_kill[c]) continue;
            for (int lv=0; lv<nverts_per_cell(); lv++)
                cells[new_nb_corners++] = vert(c, lv);
            ++new_nb_cells;
        }
        cells.resize(new_nb_corners);
    }

    void Volume::delete_vertices(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nverts());
        std::vector<bool> cells_to_kill(ncells(), false);

        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(*this, v2c, c2c);

        for (int v=0; v<nverts(); v++) {
            if (!to_kill[v]) continue;
            int cir = v2c[v];
            if (cir<0) continue; // isolated vertex
            do {
                cells_to_kill[cir/nverts_per_cell()] = true;
                cir = c2c[cir];
            } while (cir != v2c[v]);
        }
        delete_cells(cells_to_kill);

        std::vector<int> old2new;
        points.delete_points(to_kill, old2new);
        for (int &v : cells)
            v = old2new[v];
    }

    void Volume::delete_isolated_vertices()  {
        std::vector<bool> to_kill(nverts(), true);
        for (int c=0; c<ncells(); c++)
            for (int lv=0; lv<nverts_per_cell(); lv++)
                to_kill[vert(c, lv)] = false;
        delete_vertices(to_kill);
    }

    int Volume::create_cells(const int n) {
        assert(n>=0);
        cells.resize(cells.size()+n*nverts_per_cell());
        resize_attrs();
        return ncells()-n;
    }

    void compute_corner_to_corner_map(const Volume &m, std::vector<int> &v2c, std::vector<int> &c2c) {
        c2c = std::vector<int>(m.ncorners(), -1);
        v2c = std::vector<int>(m.nverts(),   -1);

        // compute vertex-to-cell-corner map
        for (int c=0; c<m.ncells(); c++)
            for (int lv=0; lv<m.nverts_per_cell(); lv++)
                v2c[m.vert(c, lv)] = m.corner(c, lv);

        // chain cell corners around vertices
        for (int c=0; c<m.ncells(); c++)
            for (int lv=0; lv<m.nverts_per_cell(); lv++) {
                int crn = m.corner(c, lv);
                int v = m.vert(c, lv);
                c2c[crn] = v2c[v];
                v2c[v] = crn;
            }
    }
}

