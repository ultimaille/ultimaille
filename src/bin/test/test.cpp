#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }

    PolyMesh pm;
    SurfaceAttributes attributes = read_geogram(argv[1], pm);
    PointAttribute<int> prand("rand", attributes, pm);
    FacetAttribute<int> fid("id", attributes, pm);
    CornerAttribute<int> cid("id", attributes, pm);

    FacetAttribute<int> nonex("nonexisting", attributes, pm);
    for (int i=0; i<pm.nfacets(); i++)
        nonex[i] = rand()%1980;

    write_geogram("read_test.geogram", pm, attributes);
    write_geogram("read_test_wo_attributes.geogram", pm);

    return 0;
}

