#include <iostream>

#include <ultimaille/all.h>

using namespace UM;

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " mesh.geogram" << std::endl;
        return 1;
    }

    PolyLine mseg;
    Triangles mtri;
    Quads mquads;
    Polygons mpoly;
    Tetrahedra mtet;
    Hexahedra mhex;

    VolumeAttributes tetattr = read_geogram(argv[1], mtet);
    CellFacetAttribute<int> rnd("rnd", tetattr, mtet);
    for (int f=0; f<mtet.nfacets(); f++)
        rnd[f] = rand()%12465;

    VolumeAttributes hexattr = read_geogram(argv[1], mhex);

    write_geogram("tet.geogram", mtet, tetattr);
    write_geogram("hex.geogram", mhex, hexattr);

    return 0;
}

