#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"


bool find_container(std::string name, std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > &attrs, std::shared_ptr<GenericAttributeContainer> &ptr) {
    for (auto &pair : attrs) {
        if (pair.first!=name) continue;
        ptr = pair.second;
        return true;
    }
    return false;
}

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }

    PolyMesh pm;
    std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > pattr, fattr, cattr;

    read_geogram(argv[1], pm, pattr, fattr, cattr);

    std::shared_ptr<GenericAttributeContainer> ptr;
    bool ret = find_container("rand", pattr, ptr);
    assert(ret);
    PointAttribute<int> pid(pm.points, ptr);

    if (auto aint = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); aint.get()!=nullptr) {
        for (int i=0; i<pm.nverts(); i++) {
            std::cerr << aint->data[i] << std::endl;
        }
    }


//    std::cerr << pattr.size() << " " << fattr.size() << " " << cattr.size() << std::endl;

    write_geogram_ascii("read_test.geogram_ascii", pm, {{"rand", pid.ptr}}, {}, {});
    return 0;
}

