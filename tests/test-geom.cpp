#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("Test triangle geom", "[geom]") {

    Triangles m;
    m.points.create_points(4);
    m.points[0] = vec3(0,0,0);
    m.points[1] = vec3(1,0,0);
    m.points[2] = vec3(0.5,0.5,0);
    m.points[3] = vec3(1.5,0.5,0);

    m.create_facets(2);
    m.vert(0, 0) = 0;
    m.vert(0, 1) = 1;
    m.vert(0, 2) = 2;
    m.vert(1, 0) = 1;
    m.vert(1, 1) = 3;
    m.vert(1, 2) = 2;

    // Check normal
    for (auto f : m.iter_facets()) {
        CHECK(f.geom<Triangle3>().normal().z > 0);
    }

    // Check bary
    auto f = m.iter_facets().begin().f;
    CHECK(std::abs((f.geom<Triangle3>().bary_verts() - vec3{0.5,0.5/3.,0}).norm()) < 1e-4);

    // Check area
    double tri_area = f.geom<Triangle3>().unsigned_area();
    // Check with formula: A = bh / 2
    double test_area = 1. /* b */ * 0.5 /* h */ * 0.5 /* divide by two */;
    CHECK(tri_area == test_area);

}

TEST_CASE("Test quad geom", "[geom]") {
    Quads m;
    m.points.create_points(6);
    m.points[0] = vec3(0,0,0);
    m.points[1] = vec3(1,0,0);
    m.points[2] = vec3(1,1,0);
    m.points[3] = vec3(0,1,0);
    m.points[4] = vec3(2,0,0);
    m.points[5] = vec3(1,2,0);

    m.create_facets(2);
    m.vert(0, 0) = 0;
    m.vert(0, 1) = 1;
    m.vert(0, 2) = 2;
    m.vert(0, 3) = 3;

    m.vert(1, 0) = 1;
    m.vert(1, 1) = 4;
    m.vert(1, 2) = 5;
    m.vert(1, 3) = 2;

    // Check normal
    for (auto f : m.iter_facets()) {
        CHECK(f.geom<Quad3>().normal().z > 0);
    }

    // Check area
    auto f = m.iter_facets().begin().f;
    double quad_area = f.geom<Quad3>().unsigned_area();
    // Considering the quad is a square of length 1 and height 1
    // We have an area equal to 1*1=1
    bool is_area_correct = quad_area >= 0.9999 && quad_area <= 1.0001;
    CHECK(is_area_correct);
}
