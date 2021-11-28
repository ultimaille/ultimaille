#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

// TODO other types of meshes

TEST_CASE("Polygons Attributes", "[Attributes]") {
    static const std::string filename = "ultimaille-test-attributes-polygons.geogram";
    Polygons m;
    *m.points.data = {{ 1 , 0 ,0}, { 0.309017, 0.951057,0}, {-0.809017, 0.587785,0}, {-0.809017,-0.587785,0}, { 0.309017,-0.951057,0}};
    m.facets = {0,1,2,3,4,0,2};
    m.offset = {0,3,7};

    PointAttribute<bool> vbool(m.points, true);
    FacetAttribute<vec3> fvec3(m, {3,14,15});
    CornerAttribute<int> cint(m, -8);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int f : range(m.nfacets()))
        fvec3[f] = m.util.bary_verts(f);

    for (int c : range(m.ncorners()))
        cint[c] = rand();

    REQUIRE(vbool.ptr->data.size() == 5);
    REQUIRE(fvec3.ptr->data.size() == 2);
    REQUIRE(cint.ptr->data.size() == 7);

    m.points.push_back({1,1,1});
    m.points.push_back({3,2,1});

    int off = m.create_facets(1, 4);

    REQUIRE(vbool.ptr->data.size() == 7);
    REQUIRE(fvec3.ptr->data.size() == 3);
    REQUIRE(cint.ptr->data.size() == 11);

    m.vert(off, 0) = 3;
    m.vert(off, 1) = 2;
    m.vert(off, 2) = m.nverts()-2;
    m.vert(off, 3) = m.nverts()-1;

    REQUIRE((fvec3[off]-vec3(3,14,15)).norm2()==0);
    for (int i : range(4)) {
        REQUIRE(cint[m.ncorners()-1-i]==-8);
    }
    REQUIRE(vbool[m.nverts()-2] == true);
    REQUIRE(vbool[m.nverts()-1] == true);

    std::vector<bool> to_kill(m.nverts(), false);
    to_kill[4] = true;
    m.delete_vertices(to_kill);

    REQUIRE(vbool.ptr->data.size() == 6);
    REQUIRE(fvec3.ptr->data.size() == 2);
    REQUIRE(cint.ptr->data.size() == 7);

    // TODO compress test

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"fvec3", fvec3.ptr}}, {{"cint", cint.ptr}}});

    Polygons m2;
    SurfaceAttributes attrs = read_by_extension(filename, m2);
    FacetAttribute<int> fint2(m2, 18);
    CornerAttribute<int> cint2("cint", attrs, m2, -122);
    PointAttribute<bool> vbool2("vbool", attrs, m2, true);

    m2.create_facets(1, 5);
    CHECK(fint2[m2.nfacets()-1]==18);
    CHECK(cint2[m2.ncorners()-1]==-122);

    int offv = m2.points.create_points(1);
    CHECK( vbool2[offv] );
}

