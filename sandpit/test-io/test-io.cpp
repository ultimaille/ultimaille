#include <iostream>

#include <ultimaille/all.h>

using namespace UM;

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " mesh.geogram" << std::endl;
        return 1;
    }
    /*
    Triangles mplg;
    SurfaceAttributes  plgattr = read_wavefront_obj(argv[1], mplg);
    PointAttribute<vec2> tex_coord("tex_coord", plgattr, mplg);
      std::cerr << tex_coord[2] << std::endl;
    */

    PolyLine mseg;
    Triangles mtri;
    Quads mqua;
    Polygons mplg;
    Tetrahedra mtet;
    Hexahedra mhex;
    Wedges mwdg;

    PolyLineAttributes segattr = read_geogram(argv[1], mseg);

    SurfaceAttributes  plgattr = read_geogram(argv[1], mplg);
    SurfaceAttributes  triattr = read_geogram(argv[1], mtri);
    SurfaceAttributes  quaattr = read_geogram(argv[1], mqua);

    VolumeAttributes   tetattr = read_geogram(argv[1], mtet);
    VolumeAttributes   hexattr = read_geogram(argv[1], mhex);
    VolumeAttributes   wdgattr = read_by_extension(argv[1], mwdg);

    CellFacetAttribute<int> rnd("rnd", tetattr, mtet);
    for (int f=0; f<mtet.nfacets(); f++)
        rnd[f] = rand()%12465;

    write_geogram("seg.geogram", mseg, segattr);
    write_by_extension("plg.geogram", mplg, plgattr);
    write_geogram("tri.geogram", mtri, triattr);
    write_geogram("qua.geogram", mqua, quaattr);
    write_geogram("tet.geogram", mtet, tetattr);
    write_geogram("hex.geogram", mhex, hexattr);
    write_geogram("wdg.geogram", mwdg, hexattr);

    return 0;
}

