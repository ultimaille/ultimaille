#include <iostream>
#include <algorithm>
#include <cassert>
#include "surface.h"
#include "attributes.h"
#include "surface_connectivity.h"

namespace UM {
    // unsigned area for a 3D triangle
    inline double unsigned_area(const vec3 &A, const vec3 &B, const vec3 &C) {
        return 0.5*cross(B-A, C-A).norm();
    }

    vec3 Surface::Util::bary_verts(const int f) const {
        vec3 ave = {0, 0, 0};
        const int nbv = m.facet_size(f);
        for (int lv=0; lv<nbv; lv++)
            ave += m.points[m.vert(f, lv)];
        return ave / static_cast<double>(nbv);
    }

    void Surface::delete_isolated_vertices()  {
        std::vector<bool> to_kill(nverts(), true);
        for (int f=0; f<nfacets(); f++)
            for (int lv=0; lv<facet_size(f); lv++)
                to_kill[vert(f, lv)] = false;
        delete_vertices(to_kill);
    }

    double Triangles::Util::unsigned_area(const int f) const {
        const vec3 &A = m.points[m.vert(f, 0)];
        const vec3 &B = m.points[m.vert(f, 1)];
        const vec3 &C = m.points[m.vert(f, 2)];
        return UM::unsigned_area(A, B, C);
    }

    void Triangles::Util::project(const int f, vec2& z0, vec2& z1, vec2& z2) const {
        const vec3 &A = m.points[m.vert(f, 0)];
        const vec3 &B = m.points[m.vert(f, 1)];
        const vec3 &C = m.points[m.vert(f, 2)];

        vec3 X = (B - A).normalized(); // construct an orthonormal 3d basis
        vec3 Z = cross(X, C - A).normalized();
        vec3 Y = cross(Z, X);

        z0 = vec2(0,0); // project the triangle to the 2d basis (X,Y)
        z1 = vec2((B - A).norm(), 0);
        z2 = vec2((C - A)*X, (C - A)*Y);
    }

    vec3 Triangles::Util::normal(const int f) const {
        const vec3 &A = m.points[m.vert(f, 0)];
        const vec3 &B = m.points[m.vert(f, 1)];
        const vec3 &C = m.points[m.vert(f, 2)];
        return cross(B-A, C-A).normalized();
    }

    double Quads::Util::unsigned_area(const int f) const {
        double a = 0;
        vec3 G = m.util.bary_verts(f);
        for (int lv=0; lv<4; lv++) {
            const vec3 &A = m.points[m.vert(f,  lv     )];
            const vec3 &B = m.points[m.vert(f, (lv+1)%4)];
            a += UM::unsigned_area(G, A, B);
        }
        return a;
    }

    vec3 Quads::Util::normal(const int f) const {
        vec3 res = {0, 0, 0};
        vec3 bary = m.util.bary_verts(f);
        for (int lv=0; lv<4; lv++)
            res += cross(
                    m.points[m.vert(f,  lv     )]-bary,
                    m.points[m.vert(f, (lv+1)%4)]-bary
                    );
        return res.normalized();
    }

    void Surface::resize_attrs() {
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->resize(nfacets());
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->resize(ncorners());
    }

    void Surface::compress_attrs(const std::vector<bool> &facets_to_kill) {
        assert(facets_to_kill.size()==(size_t)nfacets());
        std::vector<int>  facets_old2new(nfacets(),  -1);
        std::vector<int> corners_old2new(ncorners(), -1);

        int new_nb_facets  = 0;
        int new_nb_corners = 0;

        for (int f=0; f<nfacets(); f++) {
            if (facets_to_kill[f]) continue;
            for (int lv=0; lv<facet_size(f); lv++)
                corners_old2new[corner(f, lv)] = new_nb_corners++;
            facets_old2new[f] = new_nb_facets++;
        }
//      std::cerr << "compressing facet attributes\n";
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->compress(facets_old2new);
//      std::cerr << "compressing corner attributes\n";
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->compress(corners_old2new);
    }

    void Surface::delete_vertices(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nverts());
        std::vector<bool> facets_to_kill(nfacets(), false);
        SurfaceConnectivity fec(*this);

        for (int v=0; v<nverts(); v++) {
            if (!to_kill[v]) continue;
            int cir = fec.v2c[v];
            if (cir<0) continue; // isolated vertex
            do {
                facets_to_kill[fec.c2f[cir]] = true;
                cir = fec.c2c[cir];
            } while (cir != fec.v2c[v]);
        }
        delete_facets(facets_to_kill);

        std::vector<int> old2new;
        points.delete_points(to_kill, old2new);
        for (int &v : facets)
            v = old2new[v];
    }

    void Surface::delete_facets(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nfacets());
        compress_attrs(to_kill);

        int new_nb_corners = 0;
        for (int f=0; f<nfacets(); f++) {
            if (to_kill[f]) continue;
            for (int lv=0; lv<facet_size(f); lv++)
                facets[new_nb_corners++] = vert(f, lv);
        }
        facets.resize(new_nb_corners);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Triangles::create_facets(const int n) {
        facets.resize(facets.size()+n*3);
        resize_attrs();
        return nfacets()-n;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Quads::create_facets(const int n) {
        facets.resize(facets.size()+n*4);
        resize_attrs();
        return nfacets()-n;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Polygons::create_facets(const int n, const int size) {
        for (int i=0; i<n*size; i++)
            facets.push_back(0);
        for (int i=0; i<n; i++)
            offset.push_back(offset.back()+size);
        resize_attrs();
        return nfacets()-n;
    }

    void Polygons::delete_facets(const std::vector<bool> &to_kill) {
        Surface::delete_facets(to_kill);
        int new_nb_facets = 0;
        for (int f=0; f<nfacets(); f++) {
            if (to_kill[f]) continue;
            offset[new_nb_facets+1] = offset[new_nb_facets] + facet_size(f);
            ++new_nb_facets;
        }
        offset.resize(new_nb_facets+1);
    }
}

