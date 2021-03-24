#include <catch2/catch.hpp>

#include <iostream>
#include <cstdint>

#include <ultimaille/all.h>

using namespace UM;

// definition of the bunny model, check the JFIG2020 rendering contest entry: https://github.com/ssloy/tinyraytracer/wiki/Part-3:-shadertoy
constexpr int BUNW = 6;
constexpr int BUNH = 6;
constexpr int BUND = 4;
constexpr uint32_t bunny_bitfield[] = { 0xc30d0418u, 0x37dff3e0u, 0x7df71e0cu, 0x004183c3u, 0x00000400u };

static bool bunny(const int i, const int j, const int k) {
    int cubeID = i+j*BUNW+k*BUNW*BUNH;
    if (cubeID<0 || cubeID>=BUNW*BUNH*BUND) return false;
    return bunny_bitfield[cubeID/32] & (1u << (cubeID&31));
}

TEST_CASE("Hexahedra", "[VolumeConnectivity]") {
    Hexahedra m;

    { // create independent voxels
        for (int i : range(BUNW)) for (int j : range(BUNH)) for (int k : range(BUND)) {
            if (!bunny(i,j,k)) continue;
            int off_c = m.create_cells(1);
            vec3 bbox[2] = {vec3(i,j,-k), vec3(i+1,j+1,-k+1)};
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
        for (int c : cell_iter(m))
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

    VolumeConnectivity vec(m);

    CellFacetAttribute<bool> boundary(m);
    int nbrd = 0;
    for (int f : range(m.nfacets())) {
        boundary[f] = (vec.adjacent[f]<0);
        nbrd += boundary[f];
    }
    REQUIRE( nbrd==136 );

    std::vector<int> cnt_oppf(vec.max_f*vec.max_h*m.ncells(), 0);
    std::vector<int> cnt_next(vec.max_f*vec.max_h*m.ncells(), 0);
    std::vector<int> cnt_prev(vec.max_f*vec.max_h*m.ncells(), 0);

    for (int c : cell_iter(m))
        for (int lf : range(6)) {
            for (int lv : range(4)) {
                int he = vec.halfedge(c, lf, lv);
                REQUIRE( vec.cell(he)==c );
                REQUIRE( vec.facet(he)==m.facet(c, lf) );
                REQUIRE( vec.cell_facet(he)==lf );
                REQUIRE( vec.facet_halfedge(he)==lv );
                REQUIRE( vec.facet_size(he)==4 );
                REQUIRE( vec.from(he) == m.facet_vert(c, lf, lv) );
                REQUIRE( vec.to(he) == m.facet_vert(c, lf, (lv+1)%4) );
                int hen = vec.next(he);
                REQUIRE( hen>=0 );
                REQUIRE( vec.prev(hen)==he );
                cnt_next[hen]++;

                int hep = vec.prev(he);
                REQUIRE( hep>=0 );
                REQUIRE( vec.next(hep)==he );
                cnt_prev[hep]++;

                int oppc = vec.opposite_c(he);
                REQUIRE((
                        (oppc>=0 && !boundary[vec.facet(he)])
                        ||
                        (oppc==-1 && boundary[vec.facet(he)])
                       ));

                int oppf = vec.opposite_f(he);
                REQUIRE(vec.cell(oppf) == c);
                cnt_oppf[oppf]++;
            }
        }

    for (int i : range(vec.max_f*vec.max_h*m.ncells())) {
        REQUIRE( cnt_oppf[i]==1 );
        REQUIRE( cnt_next[i]==1 );
        REQUIRE( cnt_prev[i]==1 );
    }

//  write_geogram("bunny.geogram", m, {{}, {}, {{"boundary", boundary.ptr}}, {}});

/*
    Quads q;
    q.points = m.points;
    q.create_facets(nbrd);
    int cnt = 0;
    for (int f : range(m.nfacets())) {
        if (!boundary[f]) continue;
        for (int lv : range(4))
            q.vert(cnt, lv) = m.facet_vert(m.cell_from_facet(f), f%m.nfacets_per_cell(), lv);
        cnt++;
    }
    write_geogram("bunny2.geogram", q, {{}, {}, {}});
*/

}

