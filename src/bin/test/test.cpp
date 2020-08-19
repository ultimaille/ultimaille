#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }
    PolyMesh pm;
    read_geogram(argv[1], pm);

    return 0;
}

