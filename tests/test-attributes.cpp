#include <catch2/catch_test_macros.hpp>

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
    *m.points.data = {{ 1 , 0 ,0}, { 0.309017, 0.951057,0}, {-0.809017, 0.587785,0}, {-0.809017,-0.587785,0}, { 0.309017,-0.951057,0}, {0,0,0}};
    m.facets = {0,1,2,3,4,0,2};
    m.offset = {0,3,7};

    PointAttribute<bool> vbool(m.points, true);
    FacetAttribute<vec3> fvec3(m, {3,14,15});
    CornerAttribute<int> cint(m, -8);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int f : range(m.nfacets())) {
        Poly3 poly = Surface::Facet(m, f);
        fvec3[f] = poly.bary_verts();
    }

    for (int c : range(m.ncorners()))
        cint[c] = rand();

    REQUIRE(vbool.ptr->data.size() == 6);
    REQUIRE(fvec3.ptr->data.size() == 2);
    REQUIRE(cint.ptr->data.size() == 7);

    m.points.push_back({1,1,1});
    m.points.push_back({3,2,1});

#if 0
    m.connect();
    int off = m.conn->create_facet({3,2,m.nverts()-2, m.nverts()-1});
#else
    int off = m.create_facets(1, 4);
    m.vert(off, 0) = 3;
    m.vert(off, 1) = 2;
    m.vert(off, 2) = m.nverts()-2;
    m.vert(off, 3) = m.nverts()-1;
#endif

    REQUIRE(vbool.ptr->data.size() == 8);
    REQUIRE(fvec3.ptr->data.size() == 3);
    REQUIRE(cint.ptr->data.size() == 11);


    REQUIRE((fvec3[off]-vec3(3,14,15)).norm2()==0);
    for (int i : range(4)) {
        REQUIRE(cint[m.ncorners()-1-i]==-8);
    }
    REQUIRE(vbool[m.nverts()-2] == true);
    REQUIRE(vbool[m.nverts()-1] == true);

    {
        PointAttribute<int> dead_attribute(m);
        dead_attribute[3] = 14;
    }

    PointAttribute<bool> to_kill(m, false);
    to_kill[5] = true;

    REQUIRE( m.points.attr.size()==3 );
    m.delete_vertices(to_kill.ptr->data);
    REQUIRE( m.points.attr.size()==2 );


    REQUIRE(vbool.ptr->data.size() == 7);
    REQUIRE(fvec3.ptr->data.size() == 3);
    REQUIRE(cint.ptr->data.size() == 11);

    // TODO compress test

//  write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"fvec3", fvec3.ptr}}, {{"cint", cint.ptr}}});

      write_by_extension(filename, m,  {{"vbool", vbool}, {"fvec3", fvec3}, {"cint", cint}} );

    Polygons m2;
    SurfaceAttributes attrs = read_by_extension(filename, m2);
    FacetAttribute<int> fint2(m2, 18);

    CornerAttribute<int> cint2(-122);
    REQUIRE( !cint2.bound() );
    bool present = cint2.bind("cint", attrs, m2);
    CHECK( present );

    PointAttribute<bool> vbool2("vbool", attrs, m2, true);

    m2.create_facets(1, 5);
    CHECK(fint2[m2.nfacets()-1]==18);
    CHECK(cint2[m2.ncorners()-1]==-122);

    int offv = m2.points.create_points(1);
    CHECK( vbool2[offv] );
}

