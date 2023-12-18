#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("Test geom", "[geom]") {

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

    for (auto f : m.iter_facets()) {
        CHECK(f.geom<Triangle3>().normal().z > 0);
    }

    auto f = m.iter_facets().begin().f;
    CHECK(std::abs((f.geom<Triangle3>().bary_verts() - vec3{0.5,0.5/3.,0}).norm()) < 1e-4);

    Quads m2;
    m2.points.create_points(6);
    m2.points[0] = vec3(0,0,0);
    m2.points[1] = vec3(1,0,0);
    m2.points[2] = vec3(1,1,0);
    m2.points[3] = vec3(0,1,0);
    m2.points[4] = vec3(2,0,0);
    m2.points[5] = vec3(1,2,0);

    m2.create_facets(2);
    m2.vert(0, 0) = 0;
    m2.vert(0, 1) = 1;
    m2.vert(0, 2) = 2;
    m2.vert(0, 3) = 3;

    m2.vert(1, 0) = 1;
    m2.vert(1, 1) = 4;
    m2.vert(1, 2) = 5;
    m2.vert(1, 3) = 2;

    for (auto f : m2.iter_facets()) {
        CHECK(f.geom<Quad3>().normal().z > 0);
    }

}
