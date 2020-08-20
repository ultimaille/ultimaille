#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"


template <typename T> bool bind_attribute(std::string name, SurfaceAttributes attributes, PolyMesh &m, PointAttribute<T> &a) {
    for (auto &pair : std::get<0>(attributes)) {
        if (pair.first!=name) continue;
        a = PointAttribute<T>(m, pair.second);
        return true;
    }
    a = PointAttribute<T>();
    return false;
}

template <typename T> bool bind_attribute(std::string name, SurfaceAttributes attributes, PolyMesh &m, FacetAttribute<T> &a) {
    for (auto &pair : std::get<0>(attributes)) {
        if (pair.first!=name) continue;
        a = FacetAttribute<T>(m, pair.second);
        return true;
    }
    a = FacetAttribute<T>();
    return false;
}

template <typename T> bool bind_attribute(std::string name, SurfaceAttributes attributes, PolyMesh &m, CornerAttribute<T> &a) {
    for (auto &pair : std::get<0>(attributes)) {
        if (pair.first!=name) continue;
        a = CornerAttribute<T>(m, pair.second);
        return true;
    }
    a = CornerAttribute<T>();
    return false;
}

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }

    PolyMesh pm;
    SurfaceAttributes attributes = read_geogram(argv[1], pm);
    PointAttribute<int>  prand;
    FacetAttribute<int>  fid;
    CornerAttribute<int> cid;
    bind_attribute("rand", attributes, pm, prand);
    bind_attribute("id",   attributes, pm, fid);
    bind_attribute("id",   attributes, pm, cid);

    write_geogram("read_test.geogram", pm, attributes);
    write_geogram("read_test_wo_attributes.geogram", pm);

    return 0;
}

