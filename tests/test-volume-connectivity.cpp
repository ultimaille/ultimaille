#include <catch2/catch.hpp>

#include <iostream>
#include <cstdint>
#include <numeric>

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

    m.delete_isolated_vertices();
    OppositeFacet conn(m);
    CellFacetAttribute<bool> boundary(m);
    int nbrd = 0;
    for (int f : range(m.nfacets())) {
        boundary[f] = conn[f]<0;
        nbrd += boundary[f];
    }
    REQUIRE( nbrd==136 );

    REQUIRE( m.heh.nhalfedges() == 24*m.ncells() );
    std::vector<int> cnt_oppf(m.heh.nhalfedges(), 0);
    std::vector<int> cnt_next(m.heh.nhalfedges(), 0);
    std::vector<int> cnt_prev(m.heh.nhalfedges(), 0);

    for (int c : cell_iter(m)) {
        for (int lf : range(6)) {
            for (int lv : range(4)) {
                int he = lv + lf*4 + c*24;
                REQUIRE( m.heh.cell(he)==c );
                REQUIRE( m.heh.facet(he)==m.facet(c, lf) );
                REQUIRE( m.heh.cell_facet(he)==lf );
                REQUIRE( m.heh.facet_halfedge(he)==lv );
                REQUIRE( m.facet_size(m.heh.cell_facet(he))==4 );
                REQUIRE( m.heh.from(he) == m.facet_vert(c, lf, lv) );
                REQUIRE( m.heh.to(he) == m.facet_vert(c, lf, (lv+1)%4) );
                int hen = m.heh.next(he);
                REQUIRE( hen>=0 );
                REQUIRE( m.heh.prev(hen)==he );
                cnt_next[hen]++;

                int hep = m.heh.prev(he);
                REQUIRE( hep>=0 );
                REQUIRE( m.heh.next(hep)==he );
                cnt_prev[hep]++;

                int oppc = conn.opposite_c(he);
                REQUIRE((
                        (oppc>=0 && !boundary[m.heh.facet(he)])
                        ||
                        (oppc==-1 && boundary[m.heh.facet(he)])
                       ));

                int oppf = m.heh.opposite_f(he);
                REQUIRE(m.heh.cell(oppf) == c);
                cnt_oppf[oppf]++;
            }
        }
    }

    for (int i : range(24*m.ncells())) {
        REQUIRE( cnt_oppf[i]==1 );
        REQUIRE( cnt_next[i]==1 );
        REQUIRE( cnt_prev[i]==1 );
    }


    ///////////////////////////////////////////////////////////////////////////

    DisjointSet ds(m.heh.nhalfedges());

    for (int he : halfedge_iter(m)) {
        int opp = conn.opposite_c(he);
        if (opp<0) continue;
        CHECK((m.heh.to(he) == m.heh.from(opp) && m.heh.from(he)==m.heh.to(opp)));
        ds.merge(he, m.heh.opposite_f(opp));
    }

    std::vector<int> id2setid;
    std::vector<int> encounter(m.heh.nhalfedges(), 0);
    for (int h1 : halfedge_iter(m)) {
        for (int h2 : halfedge_around_edge_iter(conn, h1)) {
            encounter[h2]++;
            CHECK(ds.root(h1) == ds.root(h2));
        }
    }

    int nsets = ds.get_sets_id(id2setid);
    std::vector<int> setsize(nsets, 0);
    for (int s : id2setid) {
        setsize[s]++;
    }

    for (int h1 : halfedge_iter(m)) {
        CHECK(setsize[id2setid[h1]] == encounter[h1]);
    }

//sum_of_elems = std::accumulate(vector.begin(), vector.end(), 0);


/*
    for (int h1 : range(m.ncells()*24)) {
         for (int h2 : vec.halfedges_around_edge(h1)) {
                std::cerr << vec.from(h1) << "-" << vec.to(h1) << " : " << vec.from(h2) << "-" << vec.to(h2) << std::endl;
         }
    }
*/

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

