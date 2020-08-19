#include <iostream>
#include <algorithm>
#include <cassert>

#include "surface.h"
#include "attributes.h"

void Surface::resize_attrs() {
    for (auto &wp : attr_facets)  if (auto spt = wp.lock())
        spt->resize(nfacets());
    for (auto &wp : attr_corners) if (auto spt = wp.lock())
        spt->resize(ncorners());
}

void Surface::compress_attrs(std::vector<bool> &facets_to_kill) {
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
    for (auto &wp : attr_facets)  if (auto spt = wp.lock())
        spt->compress(facets_old2new);
    for (auto &wp : attr_corners) if (auto spt = wp.lock())
        spt->compress(corners_old2new);
}

int Surface::nverts() const {
    return points.size();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

int TriMesh::create_facets(const int n) {
    for (int i=0; i<n*3; i++)
        facets.push_back(0);
    resize_attrs();
    return nfacets()-n;
}

void TriMesh::delete_facets(std::vector<bool> &to_kill) {
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

int TriMesh::nfacets() const {
    assert(0==facets.size()%3);
    return facets.size()/3;
}

int TriMesh::ncorners() const {
    assert(0==facets.size()%3);
    return facets.size();
}

int TriMesh::facet_size(const int) const {
    return 3;
}

int TriMesh::facet_corner(const int fi, const int ci) const {
    assert(ci>=0 && ci<3 && fi>=0 && fi<nfacets());
    return fi*3 + ci;
}

int TriMesh::vert(const int fi, const int lv) const {
    assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
    return facets[fi*3 + lv];
}

int &TriMesh::vert(const int fi, const int lv) {
    assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
    return facets[fi*3 + lv];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

PolyMesh::PolyMesh() : Surface(), offset(1, 0) {}

int PolyMesh::create_facets(const int n, const int size) {
    for (int i=0; i<n*size; i++)
        facets.push_back(0);
    for (int i=0; i<n; i++)
        offset.push_back(offset.back()+size);
    resize_attrs();
    return nfacets()-n;
}

void PolyMesh::delete_facets(std::vector<bool> &to_kill) {
    assert(to_kill.size()==(size_t)nfacets());
    compress_attrs(to_kill);

    int new_nb_facets  = 0;
    int new_nb_corners = 0;
    for (int f=0; f<nfacets(); f++) {
        if (to_kill[f]) continue;
        for (int lv=0; lv<facet_size(f); lv++)
            facets[new_nb_corners++] = vert(f, lv);
        offset[++new_nb_facets] = new_nb_corners;
    }
    offset.resize(new_nb_facets+1);
    facets.resize(new_nb_corners);
}

int PolyMesh::nfacets() const {
    return static_cast<int>(offset.size())-1;
}

int PolyMesh::ncorners() const {
    return offset.back();
}

int PolyMesh::facet_size(const int fi) const {
    assert(fi>=0 && fi<nfacets());
    return offset[fi+1]-offset[fi];
}

int PolyMesh::facet_corner(const int fi, const int ci) const {
    return offset[fi]+ci;
}

int PolyMesh::vert(const int fi, const int lv) const {
    assert(fi>=0 && fi<nfacets());
    int n = facet_size(fi);
    assert(lv>=0 && lv<n);
    return facets[offset[fi]+lv];
}

int &PolyMesh::vert(const int fi, const int lv) {
    assert(fi>=0 && fi<nfacets());
    int n = facet_size(fi);
    assert(lv>=0 && lv<n);
    return facets[offset[fi]+lv];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

MeshConnectivity::MeshConnectivity(const Surface &p_m) : m(p_m) {
    int nbc = m.ncorners();
    c2f.resize(nbc, -1);
    c2c.resize(nbc, -1);
    v2c.resize(m.nverts(), -1);

    for (int f=0; f<m.nfacets(); f++) {
        for (int fc=0; fc<m.facet_size(f); fc++) {
            int c = m.facet_corner(f, fc);
            int v = m.vert(f, fc);
            c2f[c] = f;
//          c2c[c] = c;
            v2c[v] = c;
        }
    }
    for (int f=0; f<m.nfacets(); f++) { // if it ain't broke, don't fix it
        for (int fc=0; fc<m.facet_size(f); fc++) {
            int c = m.facet_corner(f, fc);
            int v = m.vert(f, fc);
            c2c[c] = v2c[v];
            v2c[v] = c;
        }
    }
}

int MeshConnectivity::from(int corner_id) const {
    int fi = c2f[corner_id];
    int lv = corner_id - m.facet_corner(fi, 0);
    return m.vert(fi, lv);
}

int MeshConnectivity::to(int corner_id) const {
    int fi = c2f[corner_id];
    int lv = corner_id - m.facet_corner(fi, 0);
    int n = m.facet_size(fi);
    return m.vert(fi, (lv+1)%n);
}

int MeshConnectivity::next(int corner_id) const {
    int fi = c2f[corner_id];
    int lv = corner_id - m.facet_corner(fi, 0);
    int n = m.facet_size(fi);
    return m.facet_corner(fi, (lv+1)%n);
}

int MeshConnectivity::prev(int corner_id) const {
    int fi = c2f[corner_id];
    int lv = corner_id - m.facet_corner(fi, 0);
    int n = m.facet_size(fi);
    return m.facet_corner(fi, (lv-1+n)%n);
}

int MeshConnectivity::opposite(int corner_id) const {
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

bool MeshConnectivity::is_border_vert(int v) {
    int cir = v2c[v];
    if (cir<0) return false;
    do {
        cir = c2c[cir];
        if (opposite(cir) == -1) return true;
    } while (cir != v2c[v]);
    return false;
}

int MeshConnectivity::next_around_vertex(int corner_id) {
    return opposite(prev(corner_id));
}

