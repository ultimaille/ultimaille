#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <vector>
#include <memory>
#include <string>
#include "geometry.h"
#include "pointset.h"
#include "attributes.h"

struct Surface { // polygonal mesh interface
    PointSet points{};
    std::vector<int> facets{};
    std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};

    Surface() = default;

    void resize_attribute_containers() {
        for (std::weak_ptr<GenericAttributeContainer> &wp : attr_facets)
            if (auto spt = wp.lock())
                spt->resize(nfacets());
    }

    void permute_attribute_containers(std::vector<int> &perm) {
        assert(false);
    }

    int nverts() const {
        return points.size();
    }

    virtual int nfacets()  const = 0;
    virtual int ncorners() const = 0;
    virtual int facet_size(const int fi) const = 0;
    virtual int facet_corner(const int fi, const int ci) const = 0;
    virtual int &vert(const int fi, const int lv) = 0;
    virtual int vert(const int fi, const int lv) const = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct TriMesh : Surface { // simplicial mesh implementation
    int create_facets(const int n, const int size);

    int nfacets()  const;
    int ncorners() const;
    int facet_size(const int fi) const;
    int facet_corner(const int fi, const int ci) const;
    int vert(const int fi, const int lv) const;
    int &vert(const int fi, const int lv);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct PolyMesh : Surface { // polygonal mesh implementation
    std::vector<int> offset;
    PolyMesh();

    int create_facets(const int n, const int size);

    int nfacets()  const;
    int ncorners() const;
    int facet_size(const int fi) const;
    int facet_corner(const int fi, const int ci) const;
    int vert(const int fi, const int lv) const;
    int &vert(const int fi, const int lv);
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

