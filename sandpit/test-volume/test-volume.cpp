#include <iostream>
#include <cstdint>

#include <ultimaille/all.h>

using namespace UM;

// definition of the bunny model, check the JFIG2020 rendering contest entry: https://github.com/ssloy/tinyraytracer/wiki/Part-3:-shadertoy
const int BUNW = 6;
const int BUNH = 6;
const int BUND = 4;
uint32_t bunny_bitfield[] = { 0xc30d0418u, 0x37dff3e0u, 0x7df71e0cu, 0x004183c3u, 0x00000400u };
bool bunny(const int i, const int j, const int k) {
    int cubeID = i+j*BUNW+k*BUNW*BUNH;
    if (cubeID<0 || cubeID>=BUNW*BUNH*BUND) return false;
    return bunny_bitfield[cubeID/32] & (1u << (cubeID&31));
}

int main() {
    Hexahedra m;

    { // create independent voxels
        for (int i : range(BUNW)) for (int j : range(BUNH)) for (int k : range(BUND)) {
            if (!bunny(i,j,k)) continue;
            int off_c = m.create_cells(1);
            vec3 bbox[2] = {vec3(i,j,-k), vec3(i,j,-k)+vec3(1,1,1)};
            for (int u : range(2)) for (int v : range(2)) for (int w : range(2)) {
                vec3 p = {bbox[u].x, bbox[v].y, bbox[w].z};
                int off_v = m.nverts();
                m.points.push_back(p);
                m.vert(off_c, u+v*2+w*4) = off_v;
            }
        }
    }

    { // colocate vertices
        std::vector<int> old2new;
        colocate(*m.points.data, old2new, 1e-3);
        for (int c : range(m.ncells()))
            for (int lv : range(8))
                m.vert(c, lv) = old2new[m.vert(c, lv)];
    }

    { // remove isolated vertices
        std::vector<bool> to_kill(m.nverts(), true);
        for (int c : range(m.ncells()))
            for (int lv : range(8))
                to_kill[m.vert(c, lv)] = false;
        m.delete_vertices(to_kill);
    }

    CellFacetAttribute<bool> boundary(m);
    { // export boundary attribute for cell facets
        VolumeConnectivity cec(m);
        for (int f : range(m.nfacets()))
            boundary[f] = (cec.adjacent[f]<0);
    }

    write_geogram("bunny.geogram", m, {{}, {}, {{"boundary", boundary.ptr}}, {}});
    return 0;
}

