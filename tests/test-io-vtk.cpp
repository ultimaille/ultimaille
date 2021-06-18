#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

static const std::string points_str =
R"(# vtk DataFile Version 4.2
written by meshio v4.3.6
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 double
0.0 0.9 0.0 0.0 0.8 0.0 1.0 0.3 0.0
)";

TEST_CASE("PointSet IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-pointset-in.vtk", "ultimaille-test-pointset-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << points_str;
    ofs.close();

    PointSet m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].size()==3 );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

static const std::string edges_str =
R"(# vtk DataFile Version 4.2
written by meshio v4.3.6
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 double
0.0 0.9 0.0 0.0 0.8 0.0 1.0 0.3 0.0
CELLS 2 6
2
0
1
2
1
2
CELL_TYPES 2
3
3
)";

TEST_CASE("Polyline IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-polyline-in.vtk", "ultimaille-test-polyline-out.vtk" };
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

static const std::string tri_str =
R"(# vtk DataFile Version 4.2
vtk output
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 double
0 0 0 1 0 0 0 1 0 

CELLS 1 4
3 0 1 2 

CELL_TYPES 1
5
)";

TEST_CASE("Triangles IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-triangles-in.vtk", "ultimaille-test-triangles-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << tri_str;
    ofs.close();

    Triangles m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==3 );
        REQUIRE( m[i].nfacets()==1 );
        CHECK( std::abs(m[i].util.unsigned_area(0)-.5)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

static const std::string quad_str =
R"(# vtk DataFile Version 4.2
vtk output
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 double
0 0 0 1 0 0 1 1 0 
0 1 0 
CELLS 1 5
4 0 1 2 3 

CELL_TYPES 1
9
)";

TEST_CASE("Quads IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-quads-in.vtk", "ultimaille-test-quads-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << quad_str;
    ofs.close();

    Quads m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].nfacets()==1 );
        CHECK( std::abs(m[i].util.unsigned_area(0)-1.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

static const std::string poly_str =
R"(# vtk DataFile Version 4.2
written by meshio v4.3.6
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 5 double
1.0 0.0 0.0 0.309017 0.951057 0.0 -0.809017 0.587785 0.0 -0.809017 -0.587785 0.0 0.309017 -0.951057 0.0
CELLS 2 9
3
0
1
2
4
3
4
0
2
CELL_TYPES 2
5
9
)";

TEST_CASE("Poly IO test", "[Polygons]") {
    static const std::string filename[2] = { "ultimaille-test-polygons-in.vtk", "ultimaille-test-polygons-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << poly_str;
    ofs.close();

    Polygons m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==5 );
        REQUIRE( m[i].nfacets()==2 );
        REQUIRE( m[i].facet_size(0)==3 );
        REQUIRE( m[i].facet_size(1)==4 );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

static const std::string tet_str =
R"(# vtk DataFile Version 4.2
vtk output
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 double
0 0 0 1 0 0 0 1 0 
0 0 1 
CELLS 1 5
4 0 1 2 3 

CELL_TYPES 1
10
)";

TEST_CASE("Tetrahedra IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-tetrahedra-in.vtk", "ultimaille-test-tetrahedra-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << tet_str;
    ofs.close();

    Tetrahedra m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].ncells()==1 );
        CHECK( std::abs(m[i].util.cell_volume(0)-1./6.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

static const std::string hex_str =
R"(# vtk DataFile Version 4.2
vtk output
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 8 double
0 0 0 1 0 0 1 1 0 
0 1 0 0 0 1 1 0 1 
1 1 1 0 1 1 
CELLS 1 9
8 0 1 2 3 4 5 6 7 

CELL_TYPES 1
12
)";

TEST_CASE("Hexahedra IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-hexahedra-in.vtk", "ultimaille-test-hexahedra-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << hex_str;
    ofs.close();

    Hexahedra m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==8 );
        REQUIRE( m[i].ncells()==1 );
        CHECK( std::abs(m[i].util.cell_volume(0)-1.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

static const std::string wedge_str =
R"(# vtk DataFile Version 4.2
vtk output
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 6 double
0 0 0 1 0 0 0 1 0 
0 0 1 1 0 1 0 1 1 

CELLS 1 7
6 0 1 2 3 4 5 

CELL_TYPES 1
13
)";

TEST_CASE("Wedges IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-wedges-in.vtk", "ultimaille-test-wedges-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << wedge_str;
    ofs.close();

    Wedges m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==6 );
        REQUIRE( m[i].ncells()==1 );
        CHECK( std::abs(m[i].util.cell_volume(0)-.5)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

static const std::string pyramid_str =
R"(# vtk DataFile Version 4.2
vtk output
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 5 double
0 0 0 1 0 0 1 1 0 
0 1 0 .5 .5 .5
CELLS 1 6
5 0 1 2 3 4 

CELL_TYPES 1
14

)";

TEST_CASE("Pyramids IO test", "[VTK]") {
    static const std::string filename[2] = { "ultimaille-test-pyramids-in.vtk", "ultimaille-test-pyramids-out.vtk" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << pyramid_str;
    ofs.close();

    Pyramids m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==5 );
        REQUIRE( m[i].ncells()==1 );
        CHECK( std::abs(m[i].util.cell_volume(0)-1./6.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}

