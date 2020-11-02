#include <iostream>
#include <algorithm>
#include <cassert>
#include "surface.h"
#include "attributes.h"

namespace UM {
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
                corners_old2new[facet_corner(f, lv)] = new_nb_corners++;
            facets_old2new[f] = new_nb_facets++;
        }
        std::cerr << "compressing facet attributes\n";
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->compress(facets_old2new);
        std::cerr << "compressing corner attributes\n";
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->compress(corners_old2new);
    }

    inline int Surface::nverts() const {
        return points.size();
    }

    inline int Surface::ncorners() const {
        return facets.size();
    }

    void Surface::delete_vertices(const std::vector<bool> &to_kill) {
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

        int new_nb_facets  = 0;
        int new_nb_corners = 0;
        for (int f=0; f<nfacets(); f++) {
            if (to_kill[f]) continue;
            for (int lv=0; lv<facet_size(f); lv++)
                facets[new_nb_corners++] = vert(f, lv);
            ++new_nb_facets;
        }
        facets.resize(new_nb_corners);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Triangles::create_facets(const int n) {
        facets.resize(facets.size()+n*3);
        resize_attrs();
        return nfacets()-n;
    }

    inline int Triangles::nfacets() const {
        assert(0==facets.size()%3);
        return facets.size()/3;
    }

    inline int Triangles::facet_size(const int) const {
        return 3;
    }

    inline int Triangles::facet_corner(const int fi, const int ci) const {
        assert(ci>=0 && ci<3 && fi>=0 && fi<nfacets());
        return fi*3 + ci;
    }

    inline int Triangles::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
        return facets[fi*3 + lv];
    }

    inline int &Triangles::vert(const int fi, const int lv) {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
        return facets[fi*3 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Quads::create_facets(const int n) {
        facets.resize(facets.size()+n*4);
        resize_attrs();
        return nfacets()-n;
    }

    inline int Quads::nfacets() const {
        assert(0==facets.size()%4);
        return facets.size()/4;
    }

    inline int Quads::facet_size(const int) const {
        return 4;
    }

    inline int Quads::facet_corner(const int fi, const int ci) const {
        assert(ci>=0 && ci<4 && fi>=0 && fi<nfacets());
        return fi*4 + ci;
    }

    inline int Quads::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<4);
        return facets[fi*4 + lv];
    }

    inline int &Quads::vert(const int fi, const int lv) {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<4);
        return facets[fi*4 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    Polygons::Polygons() : Surface(), offset(1, 0) {}

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

    void Polygons::extract_triangles(Triangles &tri) {
        tri = Triangles();
        tri.points = points;
        int ntri = 0;
        for (int f=0; f<nfacets(); f++)
            ntri += (3==facet_size(f));
        tri.create_facets(ntri);
        int cnt = 0;
        for (int f=0; f<nfacets(); f++) {
            if (3!=facet_size(f)) continue;
            for (int v=0; v<3; v++)
                tri.vert(cnt, v) = vert(f, v);
            ++cnt;
        }
    }

    void Polygons::extract_quads(Quads &quads) {
        quads = Quads();
        quads.points = points;
        int nquads = 0;
        for (int f=0; f<nfacets(); f++)
            nquads += (4==facet_size(f));
        quads.create_facets(nquads);
        int cnt = 0;
        for (int f=0; f<nfacets(); f++) {
            if (4!=facet_size(f)) continue;
            for (int v=0; v<4; v++)
                quads.vert(cnt, v) = vert(f, v);
            ++cnt;
        }
    }

    inline int Polygons::nfacets() const {
        return static_cast<int>(offset.size())-1;
    }

    inline int Polygons::facet_size(const int fi) const {
        assert(fi>=0 && fi<nfacets());
        return offset[fi+1]-offset[fi];
    }

    inline int Polygons::facet_corner(const int fi, const int ci) const {
        assert(fi>=0 && fi<nfacets());
        return offset[fi]+ci;
    }

    inline int Polygons::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets());
        int n = facet_size(fi);
        assert(lv>=0 && lv<n);
        return facets[offset[fi]+lv];
    }

    inline int &Polygons::vert(const int fi, const int lv) {
        assert(fi>=0 && fi<nfacets());
        int n = facet_size(fi);
        assert(lv>=0 && lv<n);
        return facets[offset[fi]+lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    SurfaceConnectivity::SurfaceConnectivity(const Surface &p_m) : m(p_m) {
        int nbc = m.ncorners();
        c2f.resize(nbc, -1);
        c2c.resize(nbc, -1);
        v2c.resize(m.nverts(), -1);

        for (int f=0; f<m.nfacets(); f++)
            for (int fc=0; fc<m.facet_size(f); fc++) {
                int c = m.facet_corner(f, fc);
                int v = m.vert(f, fc);
                c2f[c] = f;
                v2c[v] = c;
            }
        for (int f=0; f<m.nfacets(); f++) // if it ain't broke, don't fix it
            for (int fc=0; fc<m.facet_size(f); fc++) {
                int c = m.facet_corner(f, fc);
                int v = m.vert(f, fc);
                c2c[c] = v2c[v];
                v2c[v] = c;
            }
    }

    inline int SurfaceConnectivity::from(const int corner_id) const {
        int fi = c2f[corner_id];
        int lv = corner_id - m.facet_corner(fi, 0);
        return m.vert(fi, lv);
    }

    inline int SurfaceConnectivity::to(const int corner_id) const {
        int fi = c2f[corner_id];
        int lv = corner_id - m.facet_corner(fi, 0);
        int n = m.facet_size(fi);
        return m.vert(fi, (lv+1)%n);
    }

    inline int SurfaceConnectivity::next(const int corner_id) const {
        int fi = c2f[corner_id];
        int lv = corner_id - m.facet_corner(fi, 0);
        int n = m.facet_size(fi);
        return m.facet_corner(fi, (lv+1)%n);
    }

    inline int SurfaceConnectivity::prev(const int corner_id) const {
        int fi = c2f[corner_id];
        int lv = corner_id - m.facet_corner(fi, 0);
        int n = m.facet_size(fi);
        return m.facet_corner(fi, (lv-1+n)%n);
    }

    int SurfaceConnectivity::opposite(const int corner_id) const {
        int cir = corner_id;
        int result = -1; // not found
        do {
            int candidate = prev(cir);
            if (from(candidate) == to(corner_id) && to(candidate) == from(corner_id)) {
                if (result == -1) result = candidate;
                else return -1; // found more than one
            }
            if (cir != corner_id && to(corner_id) == to(cir))
                return -1; // the edge is non manifold
            cir = c2c[cir];
        } while (cir != corner_id);
        return result;
    }

    bool SurfaceConnectivity::is_border_vert(const int v) const {
        int cir = v2c[v];
        if (cir<0) return false;
        do {
            if (opposite(cir) == -1) return true;
            cir = c2c[cir];
        } while (cir != v2c[v]);
        return false;
    }

    int SurfaceConnectivity::next_around_vertex(const int corner_id) const {
        return opposite(prev(corner_id));
    }
}

