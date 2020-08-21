#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <vector>
#include <memory>
#include "geometry.h"
#include "pointset.h"

struct GenericAttributeContainer;

struct Surface { // polygonal mesh interface
    PointSet points{};
    std::vector<int> facets{};
    std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};
    std::vector<std::weak_ptr<GenericAttributeContainer> > attr_corners{};

    Surface() = default;

    void delete_vertices(const std::vector<bool> &to_kill);
    virtual void delete_facets(const std::vector<bool> &to_kill);
    void resize_attrs();
    void compress_attrs(const std::vector<bool> &facets_to_kill);

    int nverts() const;
    int ncorners() const;

    virtual int nfacets() const = 0;
    virtual int facet_size(const int fi) const = 0;
    virtual int facet_corner(const int fi, const int ci) const = 0;
    virtual int  vert(const int fi, const int lv) const = 0;
    virtual int &vert(const int fi, const int lv)       = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Triangles : Surface { // simplicial mesh implementation
    int create_facets(const int n);

    int nfacets()  const;
    int facet_size(const int) const;
    int facet_corner(const int fi, const int ci) const;
    int  vert(const int fi, const int lv) const;
    int &vert(const int fi, const int lv);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Quads : Surface { // quad mesh implementation
    int create_facets(const int n);

    int nfacets()  const;
    int facet_size(const int) const;
    int facet_corner(const int fi, const int ci) const;
    int  vert(const int fi, const int lv) const;
    int &vert(const int fi, const int lv);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Polygons : Surface { // polygonal mesh implementation
    std::vector<int> offset;
    Polygons();

    int create_facets(const int n, const int size);
    virtual void delete_facets(const std::vector<bool> &to_kill);
    void extract_triangles(Triangles &tri);
    void extract_quads(Quads &quads);

    int nfacets()  const;
    int facet_size(const int fi) const;
    int facet_corner(const int fi, const int ci) const;
    int  vert(const int fi, const int lv) const;
    int &vert(const int fi, const int lv);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SurfaceConnectivity { // half-edge-like connectivity interface
    SurfaceConnectivity(const Surface &p_m);

    int from(const int corner_id) const;
    int   to(const int corner_id) const;
    int prev(const int corner_id) const;
    int next(const int corner_id) const;
    int opposite(const int corner_id) const;
    bool is_border_vert(const int v) const;
    int next_around_vertex(const int corner_id) const;

    const Surface &m;
    std::vector<int> v2c; // vertex to corner
    std::vector<int> c2f; // corner to facet
    std::vector<int> c2c; // corner to next corner sharing the same vertex (unordered)
};

#endif //__SURFACE_H__

