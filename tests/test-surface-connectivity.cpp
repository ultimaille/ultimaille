#include <catch2/catch_test_macros.hpp>

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

TEST_CASE("Quads", "[SurfaceConnectivity]") {
    Quads q;
    {
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
            for (int c : m.iter_cells())
                for (int lv : range(8))
                    m.vert(c, lv) = old2new[m.vert(c, lv)];
        }

        m.delete_isolated_vertices();

        OppositeFacet conn(m);
        CellFacetAttribute<bool> boundary(m);
        int nbrd = 0;
        for (int f : range(m.nfacets())) {
            boundary[f] = conn[f]<0;
            nbrd += boundary[f];
        }

        q.points = m.points;
        q.create_facets(nbrd);
        int cnt = 0;
        for (int f : range(m.nfacets())) {
            if (!boundary[f]) continue;
            for (int lv : range(4))
                q.vert(cnt, lv) = m.facet_vert(m.cell_from_facet(f), f%m.nfacets_per_cell(), lv);
            cnt++;
        }


    }

    REQUIRE( q.nfacets()==136 );
    std::vector<bool> to_kill(136, false);
    to_kill[rand()%136] = true;
    q.delete_facets(to_kill);
    REQUIRE( q.nfacets()==135 );

//    using Vertex = typename Surface::Vertex;
    using Halfedge = typename Surface::Halfedge;
    using Facet = typename Surface::Facet;


    std::vector<int> cnt_opp(q.ncorners(), 0);
    std::vector<int> cnt_next(q.ncorners(), 0);
    std::vector<int> cnt_prev(q.ncorners(), 0);
    std::vector<int> cnt_facet(q.nfacets(), 0);

    CornerAttribute<int> val(q);
    q.connect();

    for (Halfedge he : q.iter_halfedges()) {
        Halfedge next = he.next();
        Halfedge prev = he.prev();
        Halfedge opp  = he.opposite();
        Facet f = he.facet();

        REQUIRE( prev.next() == he );
        REQUIRE( next.prev() == he );
        if (opp>=0) {
            REQUIRE(opp.opposite()==he);
            cnt_opp[opp]++;
        }

        cnt_next[next]++;
        cnt_prev[prev]++;
        cnt_facet[f]++;

        for (Halfedge cur : he.iter_sector_halfedges()) {
            (void)cur;
            val[he]++;
        }
        REQUIRE((val[he]>=1 && val[he]<=5));

    }

    int brd = 0;
    for (Halfedge he : q.iter_halfedges()) {
        REQUIRE(cnt_next[he]==1);
        REQUIRE(cnt_prev[he]==1);
        REQUIRE(cnt_opp[he]<=1);
        if (cnt_opp[he]==0) {
            REQUIRE(he.from().on_boundary());
            REQUIRE(he.to().on_boundary());
            REQUIRE(he.on_boundary());
            brd++;
        }
    }
    REQUIRE(brd==4);

//    for (auto f : q.iter_facets()) for (const auto &h : { f.halfedge(0),f.halfedge(1) }) std::cerr << h << std::endl;

    for (Facet f : q.iter_facets()) {
        REQUIRE(cnt_facet[f]==4);
    }

//  write_geogram("bunny2.geogram", q, {{}, {}, {{"val", val.ptr}}});
}

