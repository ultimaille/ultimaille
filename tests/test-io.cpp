#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

static const std::string edges_str =
R"(MeshVersionFormatted 2

Dimension 3

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
    static const std::string filename = "ultimaille-test-polyline.mesh";
    std::ofstream ofs(filename, std::ios::binary);
    ofs << edges_str;
    ofs.close();

    PolyLine m;
    read_by_extension(filename, m);

    REQUIRE( m.nverts()==3 );
    REQUIRE( m.nsegments()==2 );
}

TEST_CASE("Medit triangles IO test", "[Triangles]") {
    static const std::string filename = "ultimaille-test-triangles.mesh";
    std::ofstream ofs(filename, std::ios::binary);
    ofs << tri_str;
    ofs.close();

    Triangles m;
    read_by_extension(filename, m);

    REQUIRE( m.nverts()==3 );
    REQUIRE( m.nfacets()==1 );
    REQUIRE( std::abs(m.util.area(0)-.5)<ftol );
}

TEST_CASE("Medit quads IO test", "[Quads]") {
    static const std::string filename = "ultimaille-test-quads.mesh";
    std::ofstream ofs(filename, std::ios::binary);
    ofs << quad_str;
    ofs.close();

    Quads m;
    read_by_extension(filename, m);

    REQUIRE( m.nverts()==4 );
    REQUIRE( m.nfacets()==1 );
    REQUIRE( std::abs(m.util.area(0)-1.)<ftol );
}

TEST_CASE("Medit poly IO test", "[Polygons]") {
    REQUIRE( false );
}

TEST_CASE("Medit tetra IO test", "[Tetrahedra]") {
    static const std::string filename = "ultimaille-test-tetrahedra.mesh";
    std::ofstream ofs(filename, std::ios::binary);
    ofs << tet_str;
    ofs.close();

    Tetrahedra m;
    read_by_extension(filename, m);

    REQUIRE( m.nverts()==4 );
    REQUIRE( m.ncells()==1 );
    REQUIRE( std::abs(m.util.cell_volume(0)-1./6.)<ftol );
}

TEST_CASE("Medit hexa IO test", "[Hexahedra]") {
    static const std::string filename = "ultimaille-test-hexahedra.mesh";
    std::ofstream ofs(filename, std::ios::binary);
    ofs << hex_str;
    ofs.close();

    Hexahedra m;
    read_by_extension(filename, m);

    REQUIRE( m.nverts()==8 );
    REQUIRE( m.ncells()==1 );
    REQUIRE( std::abs(m.util.cell_volume(0)-1.)<ftol );
}

TEST_CASE("Medit wedges IO test", "[Wedges]") {
    static const std::string filename = "ultimaille-test-wedges.mesh";
    std::ofstream ofs(filename, std::ios::binary);
    ofs << wedge_str;
    ofs.close();

    Wedges m;
    read_by_extension(filename, m);

    REQUIRE( m.nverts()==6 );
    REQUIRE( m.ncells()==1 );
    REQUIRE( std::abs(m.util.cell_volume(0)-.5)<ftol );
}

TEST_CASE("Medit pyramids IO test", "[Pyramids]") {
    static const std::string filename = "ultimaille-test-pyramids.mesh";
    std::ofstream ofs(filename, std::ios::binary);
    ofs << pyramid_str;
    ofs.close();

    REQUIRE( false );
    /*
    Pyramids m;
    read_by_extension(filename, m);

    REQUIRE( m.nverts()==5 );
    REQUIRE( m.ncells()==1 );
    REQUIRE( std::abs(m.util.cell_volume(0)-1.)<ftol ); // TODO: COMPUTE THAT
    */
}

