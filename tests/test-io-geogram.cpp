#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

TEST_CASE("GEOGRAM tetrahedra + VolumeAttributes IO test", "[Tetrahedra]") {
    static const std::string filename = "ultimaille-test-tetrahedra.geogram";
    Tetrahedra m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,1,1}};
    m.cells = {0,1,2,3,4,3,2,1};

    PointAttribute<bool> vbool(m.points);
    CellAttribute<double> cdouble(m);
    CellFacetAttribute<vec3> cfvec3(m);
    CellCornerAttribute<int> ccint(m);

    for (int v : range(5))
        vbool[v] = rand()&1;

    for (int c : range(2))
        cdouble[c] = m.util.cell_volume(c);

    for (int cf : range(8))
        cfvec3[cf] = m.util.bary_facet(cf/4, cf%4);

    for (int cc : range(8))
        ccint[cc] = rand();

    write_geogram(filename, m, {{{"vbool", vbool.ptr}}, {{"cdouble", cdouble.ptr}}, {{"cfvec3", cfvec3.ptr}}, {{"ccint", ccint.ptr}}});

    Tetrahedra m2;
    VolumeAttributes attrs = read_geogram(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.ncells()==m2.ncells() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    CellAttribute<double> cdouble2("cdouble", attrs, m2);
    for (int c=0; c<m.ncells(); c++) {
        for (int lv=0; lv<m.nverts_per_cell(); lv++) {
            REQUIRE( m.vert(c, lv)==m2.vert(c, lv) );
        }
        REQUIRE( std::fabs(cdouble[c]-cdouble2[c])<ftol );
    }

    CellFacetAttribute<vec3> cfvec32("cfvec3", attrs, m2);
    for (int cf=0; cf<m.nfacets(); cf++) {
        REQUIRE( (cfvec3[cf]-cfvec32[cf]).norm()<ftol );
    }

    CellCornerAttribute<int> ccint2("ccint", attrs, m2);
    for (int cc=0; cc<m.ncorners(); cc++) {
        REQUIRE( ccint[cc]==ccint2[cc] );
    }
}

/*
static const std::string edges_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
3
0. 0.9 0. 0
0. 0.8 0. 0
1. 0.3 0. 0

Edges
2
1 2 0
2 3 0

End)";

static const std::string tri_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
3
0. 0. 0. 1
1. 0. 0. 1
0. 1. 0. 1

Triangles
1
1 2 3 1

End)";

static const std::string quad_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
4
0. 0. 0. 1
1. 0. 0. 1
1. 1. 0. 1
0. 1. 0. 1

Quadrilaterals
1
1 2 3 4  1

End
)";

static const std::string tet_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
4
0. 0. 0. 1
1. 0. 0. 1
0. 1. 0. 1
0. 1. 1. 1

Tetrahedra
1
1 2 3 4 1

End)";

static const std::string hex_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
8
0. 0. 0. 1
1. 0. 0. 1
1. 1. 0. 1
0. 1. 0. 1
0. 0. 1. 1
1. 0. 1. 1
1. 1. 1. 1
0. 1. 1. 1

Hexahedra
1
1 2 3 4 5 6 7 8 1

End)";

static const std::string wedge_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
6
0. 0. 0. 1
1. 0. 0. 1
0. 1. 0. 1
0. 0. 1. 1
1. 0. 1. 1
0. 1. 1. 1

Prisms
1
1 2 3 4 5 6 1

End)";

static const std::string pyramid_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
5
0.  0.  0.  1
1.  0.  0.  1
1.  1.  0.  1
0.  1.  0.  1
0.5 0.5 0.5 1

Pyramids
1
1 2 3 4 5 1

End)";

TEST_CASE("Medit polyline IO test", "[Polyline]") {
    static const std::string filename[2] = { "ultimaille-test-polyline-in.mesh", "ultimaille-test-polyline-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << edges_str;
    ofs.close();

    PolyLine m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==3 );
        REQUIRE( m[i].nsegments()==2 );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

TEST_CASE("Medit triangles IO test", "[Triangles]") {
    static const std::string filename[2] = { "ultimaille-test-triangles-in.mesh", "ultimaille-test-triangles-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << tri_str;
    ofs.close();

    Triangles m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==3 );
        REQUIRE( m[i].nfacets()==1 );
        REQUIRE( std::abs(m[i].util.area(0)-.5)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

TEST_CASE("Medit quads IO test", "[Quads]") {
    static const std::string filename[2] = { "ultimaille-test-quads-in.mesh", "ultimaille-test-quads-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << quad_str;
    ofs.close();

    Quads m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].nfacets()==1 );
        REQUIRE( std::abs(m[i].util.area(0)-1.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

TEST_CASE("Medit poly IO test", "[Polygons]") {
    REQUIRE( false );
}

TEST_CASE("Medit tetra IO test", "[Tetrahedra]") {
    static const std::string filename[2] = { "ultimaille-test-tetrahedra-in.mesh", "ultimaille-test-tetrahedra-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << tet_str;
    ofs.close();

    Tetrahedra m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].ncells()==1 );
        REQUIRE( std::abs(m[i].util.cell_volume(0)-1./6.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

TEST_CASE("Medit hexa IO test", "[Hexahedra]") {
    static const std::string filename[2] = { "ultimaille-test-hexahedra-in.mesh", "ultimaille-test-hexahedra-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << hex_str;
    ofs.close();

    Hexahedra m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==8 );
        REQUIRE( m[i].ncells()==1 );
        REQUIRE( std::abs(m[i].util.cell_volume(0)-1.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}


TEST_CASE("Medit wedges IO test", "[Wedges]") {
    static const std::string filename[2] = { "ultimaille-test-wedges-in.mesh", "ultimaille-test-wedges-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << wedge_str;
    ofs.close();

    Wedges m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==6 );
        REQUIRE( m[i].ncells()==1 );
        REQUIRE( std::abs(m[i].util.cell_volume(0)-.5)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

TEST_CASE("Medit pyramids IO test", "[Pyramids]") {
    static const std::string filename[2] = { "ultimaille-test-pyramids-in.mesh", "ultimaille-test-pyramids-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << pyramid_str;
    ofs.close();

    Pyramids m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==5 );
        REQUIRE( m[i].ncells()==1 );
        REQUIRE( std::abs(m[i].util.cell_volume(0)-1./6.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}


*/
