#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"
#include "ultimaille/range.h"
/*
inline auto facets(Surface &m) {
    return range(m.nfacets());
}

inline auto corners(Surface &m) {
    return range(m.ncorners());
}

struct facet_vertices {
    Surface &m_;
    const int facet_;
    facet_vertices(Surface &m, const int facet) : m_(m), facet_(facet) {}

    struct iterator {
        void operator++() { ++value_; }
        int &operator*() const { return m_.vert(facet_, value_); }
        bool operator!=(const iterator& rhs) const { return value_ != rhs.value_; }
        Surface &m_;
        const int facet_;
        int value_;
    };

    iterator begin() const { return iterator{m_, facet_, 0};    }
    iterator end()   const { return iterator{m_, facet_, m_.facet_size(facet_)}; }
};
*/

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }

    Polygons pm;
    SurfaceAttributes attributes = read_geogram(argv[1], pm);

/*

    Triangles tri;
    pm.extract_triangles(tri);
    assert(tri.nfacets() == pm.nfacets());

    Quads quads;
    pm.extract_quads(quads);
//  assert(quads.nfacets() == pm.nfacets());
//  std::cerr << quads.points.use_count() << std::endl;

//  PointAttribute<int> prand("rand", attributes, tri);
//  FacetAttribute<int> fid("id", attributes, tri);
//  CornerAttribute<int> cid("id", attributes, tri);

    FacetAttribute<int> nonex("nonexisting", attributes, tri);
    for (int i=0; i<pm.nfacets(); i++)
        nonex[i] = rand()%1980;


*/
//  for (int f : facets(pm))
//      for (int &v : facet_vertices(pm, f))
//          v = rand()%999;

/*
    for (int f : facets(pm))
        for (auto [lv, v] : enumerate(facet_vertices(pm, f))) {
            std::cerr << lv << " " << v << std::endl;
        }
*/
    write_geogram("read_test.geogram", pm, attributes);
    return 0;
}

