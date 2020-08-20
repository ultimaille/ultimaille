#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"


template <typename A,typename M> bool bind_attribute(std::string name, std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > &attrs, M& m, A &a) {
    for (auto &pair : attrs) {
        if (pair.first!=name) continue;
        a = A(m, pair.second);
        return true;
    }
    a = A();
    return false;
}

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }

    PolyMesh pm;
    auto [pattr, fattr, cattr] = read_geogram(argv[1], pm);
    PointAttribute<int>  prand;
    FacetAttribute<int>  fid;
    CornerAttribute<int> cid;
    bind_attribute("rand", pattr, pm, prand);
    bind_attribute("id",   fattr, pm, fid);
    bind_attribute("id",   cattr, pm, cid);

    write_geogram("read_test.geogram", pm, {{"rand", prand.ptr}}, {{"id", fid.ptr}}, {{"id", cid.ptr}});

    return 0;
}

