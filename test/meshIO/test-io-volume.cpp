#include <iostream>

#include <ultimaille/mesh_io.h>
#include <ultimaille/volume.h>
#include <ultimaille/attributes.h>

using namespace UM;

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " mesh.geogram" << std::endl;
        return 1;
    }

    Tetrahedra mtet;
    Hexahedra mhex;
    VolumeAttributes tetattr = read_geogram(argv[1], mtet);
    VolumeAttributes hexattr = read_geogram(argv[1], mhex);

    write_geogram("tet.geogram", mtet, tetattr);
    write_geogram("hex.geogram", mhex, hexattr);

    return 0;
}

