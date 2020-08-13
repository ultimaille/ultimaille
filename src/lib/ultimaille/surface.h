#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <vector>
#include <list>
#include <string>
#include "geometry.h"

struct GenericAttribute {
    GenericAttribute() = default;
    std::string name{};

    virtual void resize(int n) = 0;
    virtual void permute(std::vector<int> &old2new) = 0;
};

template <typename T> struct Attribute : GenericAttribute {
    void resize(int n) {
        data.resize(n);
    }
    std::vector<T> data{};
};

struct Surface { // polygonal mesh interface
    std::vector<vec3> verts{};
    std::vector<int> facets{};
    std::list<GenericAttribute *> attr_vertices{};

    Surface() = default;

    int nverts() const;

    virtual int nfacets()  const = 0;
    virtual int ncorners() const = 0;
    virtual int facet_size(const int fi) const = 0;
    virtual int facet_corner(const int fi, const int ci) const = 0;
    virtual int vert(const int fi, const int lv) const = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct TriMesh : Surface { // simplicial mesh implementation
    int nfacets()  const;
    int ncorners() const;
    int facet_size(const int fi) const;
    int facet_corner(const int fi, const int ci) const;
    int vert(const int fi, const int lv) const;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct PolyMesh : Surface { // polygonal mesh implementation
    std::vector<int> offset;

    PolyMesh();

    int nfacets()  const;
    int ncorners() const;
    int facet_size(const int fi) const;
    int facet_corner(const int fi, const int ci) const;
    int vert(const int fi, const int lv) const;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct MeshConnectivity { // half-edge-like connectivity interface
    MeshConnectivity(const Surface &p_m);

    int from(int corner_id) const;
    int   to(int corner_id) const;
    int prev(int corner_id) const;
    int next(int corner_id) const;
    int opposite(int corner_id) const;
    bool is_border_vert(int v);
    int next_around_vertex(int corner_id);

    const Surface &m;
    std::vector<int> v2c; // vertex to corner
    std::vector<int> c2f; // corner to facet
    std::vector<int> c2c; // corner to next corner sharing the same vertex (unordered)
};

#endif //__SURFACE_H__

