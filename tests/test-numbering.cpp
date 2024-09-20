#include <catch2/catch_test_macros.hpp>

#include <iostream>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

/*
TEST_CASE("Triangles numbering convention test", "[numbering convention]") {
    CHECK( false );
}

TEST_CASE("Quads numbering convention test", "[numbering convention]") {
    CHECK( false );
}

TEST_CASE("Polygons numbering convention test", "[numbering convention]") {
    CHECK( false );
}
*/

TEST_CASE("Tetrahedra numbering convention test", "[numbering convention]") {
    Tetrahedra m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,1,1}};
    m.cells = {0,1,2,3,4,3,2,1};

    // positive volume for the right-hand orientation
    Volume::Cell c0(m, 0);
    CHECK( std::abs(Tetrahedron(c0).volume()-1./6.)<ftol );
    CHECK( std::abs(Tetrahedron(Volume::Cell(m, 1)).volume()-1./3.)<ftol );

    // normals pointing outside + facet i is opposite to vertex i
    const vec3 ref_nrm[] = {vec3{1,1,1}.normalized(), {-1,0,0}, {0,-1,0}, {0,0,-1}};
    for (int f : range(4)) {
        vec3 n = Triangle3(c0.facet(f)).normal();
        CHECK( (n-ref_nrm[f]).norm()<ftol );
    }

    // smallest vertex starts each cell facet
    for (int f : range(4)) {
        for (int v:range(2)) {
            CHECK( m.facet_vert(0, f, 0)<m.facet_vert(0, f, v+1) );
        }
    }
}

TEST_CASE("Hexahedra numbering convention test", "[numbering convention]") {
    Hexahedra m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
    m.cells = {0,1,2,3,4,5,6,7};

    Volume::Cell c0(m, 0);
    // positive volume for the right-hand orientation
    CHECK( std::abs(Hexahedron(c0).volume()-1.)<ftol );

    // normals pointing outside
    const vec3 ref_nrm[] = {{-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};
    for (int f : range(6)) {
        vec3 n = Quad3(c0.facet(f)).normal();
        REQUIRE( (n-ref_nrm[f]).norm()<ftol );
    }

    // smallest vertex starts each cell facet
    for (int f : range(6)) {
        for (int v:range(3)) {
            CHECK( m.facet_vert(0, f, 0)<m.facet_vert(0, f, v+1) );
        }
    }
}

TEST_CASE("Wedges numbering convention test", "[numbering convention]") {
    Wedges m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {0,1,1}};
    m.cells = {0,1,2,3,4,5};

    Volume::Cell c0(m, 0);

    // positive volume for the right-hand orientation
    CHECK( std::abs(Wedge(c0).volume()-.5)<ftol );

    // normals pointing outside
    const vec3 ref_nrm[] = {{0,0,-1}, {0,0,1}, {0,-1,0}, {-1,0,0}, {1/std::sqrt(2.),1/std::sqrt(2.),0}};
    for (int f : range(5)) {
        Volume::Facet fa = c0.facet(f);
        vec3 n = Poly3(fa).normal();
        REQUIRE( (n-ref_nrm[f]).norm()<ftol );
    }

    // smallest vertex starts each cell facet
    for (int f : range(5)) {
        for (int v:range(2)) {
            CHECK( m.facet_vert(0, f, 0)<m.facet_vert(0, f, v+1) );
        }
    }
}

TEST_CASE("Pyramids numbering convention test", "[numbering convention]") {
    Pyramids m;
    *m.points.data = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {.5,.5,.5}};
    m.cells = {0,1,2,3,4};

    Volume::Cell c0(m, 0);

    // positive volume for the right-hand orientation
    CHECK( std::abs(Pyramid(c0).volume()-1./6.)<ftol );

    // normals pointing outside
    const vec3 ref_nrm[] = {{0,0,-1}, vec3{0,-1,1}.normalized(), vec3{-1,0,1}.normalized(), vec3{0,1,1}.normalized(), vec3{1,0,1}.normalized()};
    for (int f : range(5)) {
        Volume::Facet fa = c0.facet(f);
        vec3 n = Poly3(fa).normal();
        REQUIRE( (n-ref_nrm[f]).norm()<ftol );
    }
    // smallest vertex starts each cell facet
    for (int f : range(5)) {
        for (int v:range(2)) {
            CHECK( m.facet_vert(0, f, 0)<m.facet_vert(0, f, v+1) );
        }
    }
}

