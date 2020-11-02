#include <iostream>
#include <cstdint>

#include <ultimaille/all.h>

using namespace UM;

// definition of the bunny model, check the rendering contest entry: https://github.com/ssloy/tinyraytracer/wiki/Part-3:-shadertoy

const int BUNW = 6;
const int BUNH = 6;
const int BUND = 4;

uint32_t bunny_bitfield[] = { 0xc30d0418u, 0x37dff3e0u, 0x7df71e0cu, 0x004183c3u, 0x00000400u };
bool bunny(const int cubeID) {
    if (cubeID<0 || cubeID>=BUNW*BUNH*BUND) return false;
    return 0u != (bunny_bitfield[cubeID/32] & (1u << (cubeID&31)));
}

int main(int argc, char** argv) {
    Hexahedra m;

    // create independent voxels
    for (int i : range(BUNW)) for (int j : range(BUNH)) for (int k : range(BUND)) {
        int cellID = i+j*BUNW+k*BUNW*BUNH;
        if (!bunny(cellID)) continue;
        int off_c = m.create_hexa(1);
        vec3 bbox[2] = {vec3(i,j,-k), vec3(i,j,-k)+vec3(1,1,1)};
        for (int u : range(2)) for (int v : range(2)) for (int w : range(2)) {
            vec3 p = {bbox[u].x, bbox[v].y, bbox[w].z};
            int off_v = m.nverts();
            m.points.push_back(p);
            m.vert(off_c, u+v*2+w*4) = off_v;
        }
    }

    // colocate vertices

    std::vector<int> old2new;
    colocate(*m.points.data, old2new, 1e-3);
    for (int c : range(m.ncells()))
        for (int lv : range(8))
            m.vert(c, lv) = old2new[m.vert(c, lv)];

    // TODO remove isolated vertices

    VolumeConnectivity fec(m);
    CellFacetAttribute<bool> boundary(m);

    for (int f : range(m.nfacets()))
        boundary[f] = (fec.adjacent[f]<0);

    write_geogram("bunny.geogram", m, {{}, {}, {{"boundary", boundary.ptr}}, {}});

    return 0;
}

