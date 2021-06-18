#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

static const std::string points_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
3
0. 0.9 0. 0
0. 0.8 0. 0
1. 0.3 0. 0

End
)";

TEST_CASE("Medit PointSet IO test", "[Medit]") {
    static const std::string filename[2] = { "ultimaille-test-pointset-in.mesh", "ultimaille-test-pointset-out.mesh" };
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

End
)";

TEST_CASE("Medit Polyline IO test", "[Medit]") {
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

End
)";

TEST_CASE("Medit Triangles IO test", "[Medit]") {
    static const std::string filename[2] = { "ultimaille-test-triangles-in.mesh", "ultimaille-test-triangles-out.mesh" };
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

TEST_CASE("Medit Quads IO test", "[Medit]") {
    static const std::string filename[2] = { "ultimaille-test-quads-in.mesh", "ultimaille-test-quads-out.mesh" };
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
R"(MeshVersionFormatted 2

Dimension
3

Vertices
5
 1         0        0 0
 0.309017  0.951057 0 0
-0.809017  0.587785 0 0
-0.809017 -0.587785 0 0
 0.309017 -0.951057 0 0

Triangles
1
1 2 3 0

Quadrilaterals
1
4 5 1 3 0

End
)";

TEST_CASE("Medit poly IO test", "[Polygons]") {
    static const std::string filename[2] = { "ultimaille-test-polygons-in.mesh", "ultimaille-test-polygons-out.mesh" };
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

End
)";

TEST_CASE("Medit Tetrahedra IO test", "[Medit]") {
    static const std::string filename[2] = { "ultimaille-test-tetrahedra-in.mesh", "ultimaille-test-tetrahedra-out.mesh" };
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

End
)";

TEST_CASE("Medit Hexahedra IO test", "[Medit]") {
    static const std::string filename[2] = { "ultimaille-test-hexahedra-in.mesh", "ultimaille-test-hexahedra-out.mesh" };
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

End
)";

TEST_CASE("Medit Wedges IO test", "[Medit]") {
    static const std::string filename[2] = { "ultimaille-test-wedges-in.mesh", "ultimaille-test-wedges-out.mesh" };
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

End
)";

TEST_CASE("Medit Pyramids IO test", "[Medit]") {
    static const std::string filename[2] = { "ultimaille-test-pyramids-in.mesh", "ultimaille-test-pyramids-out.mesh" };
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

