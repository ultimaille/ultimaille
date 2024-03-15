#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch_test_macros.hpp>

#include <iostream>
#include <cstdint>
#include <numeric>
#include <iterator>

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
    
    m.connect();

    CellFacetAttribute<bool> boundary(m);
    int nbrd = 0;
    for (auto f : m.iter_facets()) {
        boundary[f] = !f.opposite().active();
        nbrd += boundary[f];
    }

    REQUIRE( nbrd==136 );

    int nhalfedges = 0;
    for (auto h : m.iter_halfedges())
        nhalfedges++;

    REQUIRE(nhalfedges == 24*m.ncells());
    

    std::vector<int> cnt_oppf(nhalfedges, 0);
    std::vector<int> cnt_next(nhalfedges, 0);
    std::vector<int> cnt_prev(nhalfedges, 0);

    for (auto c : m.iter_cells()) {
        for (auto f: c.iter_facets()) {
            for (auto h : f.iter_halfedges()) {
                
                int lf = f.id_in_cell();
                int lv = h.id_in_facet();

                REQUIRE(h.cell() == c);
                REQUIRE(h.facet() == m.facet(c, lf));
                REQUIRE(h.facet().id_in_cell() == lf);
                REQUIRE(h.id_in_facet() == lv);
                REQUIRE(m.facet_size(f) == 4);
                REQUIRE(h.from() == m.facet_vert(c, lf, lv));
                REQUIRE(h.to() == m.facet_vert(c, lf, (lv+1)%4));

                auto hen = h.next();
                REQUIRE(hen >= 0);
                REQUIRE(hen.prev() == h);
                cnt_next[hen]++;

                auto hep = h.prev();
                REQUIRE(hep >= 0);
                REQUIRE(hep.next() == h);
                cnt_prev[hep]++;

                int oppc = h.opposite_c();

                REQUIRE((
                            (oppc>=0 && !boundary[h.facet()])
                            ||
                            (oppc==-1 && boundary[h.facet()])
                        ));

                auto oppf = h.opposite_f();
                REQUIRE(oppf.cell() == c);
                cnt_oppf[oppf]++;
            }
        }
    }

    for (int i : range(24*m.ncells())) {
        REQUIRE(cnt_oppf[i] == 1);
        REQUIRE(cnt_next[i] == 1);
        REQUIRE(cnt_prev[i] == 1);
    }

    ///////////////////////////////////////////////////////////////////////////

    DisjointSet ds(nhalfedges);

    for (auto he : m.iter_halfedges()) {
        auto opp = he.opposite_c();
        if (opp<0) continue;
        CHECK((he.to() == opp.from() && he.from() == opp.to()));
        ds.merge(he, opp.opposite_f());
    }

    std::vector<int> id2setid;
    std::vector<int> encounter(nhalfedges, 0);
    for (auto h1 : m.iter_halfedges()) {
        for (auto h2 : h1.iter_CCW_around_edge()) {
            encounter[h2]++;
            CHECK(ds.root(h1) == ds.root(h2));
        }
    }

    int nsets = ds.get_sets_id(id2setid);
    std::vector<int> setsize(nsets, 0);
    for (int s : id2setid) {
        setsize[s]++;
    }

    for (auto h1 : m.iter_halfedges()) {
        CHECK(setsize[id2setid[h1]] == encounter[h1]);
    }

    // TODO see
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

#define FOR(i, n) for(int i = 0; i < (int) n; i++)

void test_volume(Volume& m) {
    m.connect();

    {
        bool facet_cell_link = true;
        for (auto c : m.iter_cells()) for (auto f : c.iter_facets())
            if (c != f.cell()) facet_cell_link = false;
        CHECK(facet_cell_link);
    }


    {
        bool opposite_facet_is_symmetric = true;
        for (auto f : m.iter_facets())
            if (f.opposite().active())
                if (f != f.opposite().opposite()) opposite_facet_is_symmetric = false;
        CHECK(opposite_facet_is_symmetric);
    }

    {
        bool opposite_c_is_symmetric = true;
        for (auto h : m.iter_halfedges())
            if (h.opposite_c().active())
                if (h != h.opposite_c().opposite_c()) opposite_c_is_symmetric = false;
        CHECK(opposite_c_is_symmetric);
    }

    {
        bool opposite_f_is_symmetric = true;
        for (auto h : m.iter_halfedges())
            if (h.opposite_f().active())
                if (h != h.opposite_f().opposite_f()) opposite_f_is_symmetric = false;
        CHECK(opposite_f_is_symmetric);
    }

    {
        bool next_prev_is_symmetric = true;
        for (auto h : m.iter_halfedges())
            if (h.prev().next() != h) next_prev_is_symmetric = false;
        for (auto h : m.iter_halfedges())
            if (h.next().prev() != h) next_prev_is_symmetric = false;
        CHECK(next_prev_is_symmetric);
    }

    {
        bool coherent_number_adjacent_cells_and_halfedges = true;
        PointAttribute<int> nb_neig_cells(m, 0);
        PointAttribute<int> nb_neig_halfedges(m, 0);

        for (auto c : m.iter_cells()) 
            for (int lv = 0; lv < c.nverts(); lv++)
                nb_neig_cells[c.vertex(lv)]++;

        for (auto h : m.iter_halfedges()) 
            nb_neig_halfedges[h.from()]++;

        for (auto v : m.iter_vertices())
            if (3 * nb_neig_cells[v] != nb_neig_halfedges[v]) 
                coherent_number_adjacent_cells_and_halfedges = false;

        CHECK(coherent_number_adjacent_cells_and_halfedges);
    }

    {
        bool halfedge_to_edge_is_coherent = true;
        bool edge_valence_is_ok = true;

        EdgeGraph edges(m);
        EdgeAttribute<int> valence(edges, 0);
        
        for (auto h : m.iter_halfedges())
            valence[edges.edge_from_halfedge(h)]++;

        for (auto e : edges.iter_edges()) {
            auto h = edges.halfedge_from_edge(e);
            int val = 0;
            
            for (auto cir : h.iter_CCW_around_edge()) {
                if (edges.edge_from_halfedge(h) != e)
                    halfedge_to_edge_is_coherent = false;
                val++;
            }

            if (val != valence[e])  edge_valence_is_ok = false;
        }

        CHECK(halfedge_to_edge_is_coherent);
        CHECK(edge_valence_is_ok);
    }

    {
        bool halfedge_to_corner_coherence = true;

        for (auto h : m.iter_halfedges())
            if (h.from_corner() != h.prev().opposite_f().from_corner()) halfedge_to_corner_coherence = false;
        
        CHECK(halfedge_to_corner_coherence);
    }
}


void split_hex(Hexahedra& hex) {

    int off_c = hex.ncells();
    int off_v = hex.points.create_points(64 * hex.ncells());
    hex.create_cells(8 * hex.ncells());

    int c = off_c;

    FOR(old_c, off_c) {

        FOR(off_i, 2) FOR(off_j, 2) FOR(off_k, 2) {
            FOR(di, 2) FOR(dj, 2) FOR(dk, 2) {

                int lv = di + 2 * dj + 4 * dk;
                int v = off_v + 8 * (c - off_c) + lv;
                hex.vert(c, lv) = v;

                double U[3] = { double(off_i + di) / 2.,double(off_j + dj) / 2.,double(off_k + dk) / 2. };
                hex.points[v] = vec3(0, 0, 0);

                FOR(old_di, 2) FOR(old_dj, 2) FOR(old_dk, 2) {
                    int old_lv = old_di + 2 * old_dj + 4 * old_dk;
                    int old_v = hex.vert(old_c, old_lv);
                    double wu = (old_di == 1) ? U[0] : (1. - U[0]);
                    double wv = (old_dj == 1) ? U[1] : (1. - U[1]);
                    double ww = (old_dk == 1) ? U[2] : (1. - U[2]);
                    hex.points[v] += wu * wv * ww * hex.points[old_v];
                }
            }

            c++;
        }
    }

    std::vector<bool> to_kill(off_c, true);
    to_kill.resize(hex.ncells(), false);
    hex.delete_cells(to_kill);

    std::vector<int> old2new;

    UM::colocate(*(hex.points.data), old2new, 1e-10);
    std::vector<bool> tokill(hex.points.size());

    FOR(v, hex.points.size()) 
        tokill[v] = (old2new[v] != v);
    
    std::vector<int> old2new2;
    hex.points.delete_points(tokill, old2new2);
    
    for (int& it : old2new) 
        it = old2new2[it];

    FOR(nc, hex.ncells()) 
        FOR(lv, 8) 
            hex.vert(nc, lv) = old2new[hex.vert(nc, lv)];
}

void hex_grid(Hexahedra& m, int n) {

    int n1 = n + 1;
    m.points.create_points(n1 * n1 * n1);

    int n2 = n1 * n1;

    FOR(i, n1 * n1 * n1) 
        m.points[i] = (1. / double(n1)) * vec3(i / n2, (i - (i / n2) * n2) / n1, i % n1);

    m.create_cells(n * n * n);

    FOR(i, n) FOR(j, n) FOR(k, n) {
        int c = i + n * j + n * n * k;

        FOR(di, 2) FOR(dj, 2) FOR(dk, 2) 
            m.vert(c, di + 2 * dj + 4 * dk) = i + di + n1 * (j + dj) + n1 * n1 * (k + dk);
    }
}

void edge_one_ring(Tetrahedra& m, int n) {

    m.points.create_points(n + 2);

    FOR(v, n) m.points[v] = {std::cos(2. * M_PI * double(v) / double(n)), std::sin(2. * M_PI * double(v) / double(n)), 0.};

    m.points[n] = vec3(0, 0, 1);
    m.points[n + 1] = vec3(0, 0, -1);
    m.create_cells(n);

    FOR(i, n) {
        m.vert(i, 0) = i;
        m.vert(i, 1) = (i + 1) % n;
        m.vert(i, 2) = n;
        m.vert(i, 3) = n + 1;
    }
}


TEST_CASE("Hex empty", "[volume][connectivity]") {
    Hexahedra hex;
    test_volume(hex);
}

TEST_CASE("Hex irregular edge", "[volume][connectivity]") {
    Hexahedra hex;
    hex.points.create_points(14);
    hex.create_cells(3);
    vec3 pts[7] = { vec3(0, 0, 0) ,vec3(.5, 0, 0) ,vec3(1, 0, 0) ,vec3(.2, .5, 0) ,vec3(.5, .4, 0) ,vec3(.8, .5, 0) ,vec3(.5, 1, 0) };
    
    FOR(i, 7) hex.points[i] = pts[i];
    FOR(lv, 7) hex.points[7 + lv] = hex.points[lv] + vec3(0, 0, .7);

    int topo[3][4] = { {0,1,3,4},{1,2,4,5},{3,4,6,5} };

    FOR(q, 3) FOR(lv, 4) hex.vert(q, lv) = topo[q][lv];
    FOR(c, 3) FOR(lv, 4) hex.vert(c, lv + 4) = hex.vert(c, lv) + 7;
    FOR(i, 3)split_hex(hex);

    test_volume(hex);
}

TEST_CASE("Hex singular vertex 1", "[volume][connectivity]") {
    Hexahedra hex;
    hex_grid(hex, 2);
    std::vector<bool> to_kill = { true,true,true,false,false,true,true,true };
    hex.delete_cells(to_kill);
    hex.delete_isolated_vertices();

    test_volume(hex);
}
TEST_CASE("Hex singular vertex 2", "[volume][connectivity]") {
    Hexahedra hex;
    hex_grid(hex, 2);
    std::vector<bool> to_kill = { false,false,false,true,true,false,false,false };
    hex.delete_cells(to_kill);
    hex.delete_isolated_vertices();
    test_volume(hex);
}
TEST_CASE("Hex singular edge", "[volume][connectivity]") {
    Hexahedra hex;
    hex_grid(hex, 2);
    std::vector<bool> to_kill = { true,true,false,false,false,false,true,true };
    hex.delete_cells(to_kill);
    hex.delete_isolated_vertices();
    test_volume(hex);
}

TEST_CASE("Tet one-ring edge", "[volume][connectivity]") {
    Tetrahedra tet;
    edge_one_ring(tet, 8);
    test_volume(tet);
}

TEST_CASE("Tet double singular edge", "[volume][connectivity]") {
    Tetrahedra tet;
    edge_one_ring(tet, 8);
    std::vector<bool> to_kill = { true, false, true, false, false, false, false, false };;
    tet.delete_cells(to_kill);
    tet.delete_isolated_vertices();
    test_volume(tet);
}

TEST_CASE("Tet triple singular edge", "[volume][connectivity]") {
    Tetrahedra tet;
    edge_one_ring(tet, 8);
    std::vector<bool> to_kill = { true, false, true, false, true, false, false, false };;
    tet.delete_cells(to_kill);
    tet.delete_isolated_vertices();
    test_volume(tet);
}

TEST_CASE("Tet singular vertex", "[volume][connectivity]") {
    Tetrahedra tet;
    tet.points.create_points(7); // no geometry here
    tet.create_cells(2);
    int verts[2][4] = { {0,1,2,3},{0,4,5,6} };
    FOR(c, 2) FOR(v, 4) tet.vert(c, v) = verts[c][v];
    test_volume(tet);
}