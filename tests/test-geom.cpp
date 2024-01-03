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
    auto tri_f = f.geom<Triangle3>();

    CHECK(std::abs((tri_f.bary_verts() - vec3{0.5,0.5/3.,0}).norm()) < 1e-4);

    // Check area comparing with formula A = bh / 2
    double tri_area = tri_f.unsigned_area();
    double test_area = 1. /* b */ * 0.5 /* h */ * 0.5 /* divide by two */;
    CHECK(std::abs(tri_area - test_area) < 1e-4);

    // Check projection
    vec2 a, b, c;
    tri_f.project(a, b, c);
    // a is near from (0,0)
    CHECK(std::abs(a.norm2()) < 1e-4);
    // b is near from the norm of b-a
    vec2 x{sqrt(1), 0};
    CHECK(std::abs((x - b).norm2()) < 1e-4);
    // c is somewhere

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

    // Get geometry of first face
    auto f = m.iter_facets().begin().f;
    auto quad_f = f.geom<Quad3>();

    INFO("quad: " << quad_f.normal());

    // Check normal
    for (auto f : m.iter_facets()) {
        auto n = quad_f.normal().z;
        CHECK(n > 0);
    }

    // Check area
    double quad_area = quad_f.unsigned_area();
    // Considering the quad is a square of length 1 and height 1
    // We have an area equal to 1*1=1
    bool is_area_correct = quad_area >= 0.9999 && quad_area <= 1.0001;
    CHECK(is_area_correct);
}

TEST_CASE("Test poly geom", "[geom]") {

    // Create a unique regular polygon facet of size 6
    const int nbv = 6;
    Polygons m;
    m.points.create_points(nbv);
    
    // Create a regular polygon
    for (int i = 0; i < nbv; i++) {
        const double t = 2*M_PI/nbv*i;
        m.points[i] = vec3(cos(t), sin(t), 0);
    }

    // Associate each vertex with corresponding point
    m.create_facets(1, nbv);
    for (int i = 0; i < nbv; i++) {
        m.vert(0, i) = i;
    }

    // // Get geometry of first face
    auto f = m.iter_facets().begin().f;
    auto poly_f = f.geom<Poly3>();

    // // Check normal
    for (auto f : m.iter_facets()) {
        CHECK(f.geom<Poly3>().normal().z > 0);
    }

    // Check bary
    vec3 bary = poly_f.bary_verts();
    CHECK(std::abs((bary - vec3{0,0,0}).norm2()) < 1e-4);

    INFO("poly computesize: " << poly_f.v.size());

    // Check area
    double area = poly_f.unsigned_area();
    INFO("polygon area is: " << area);
    // Considering the polygon is regular we have the formula: A = r²n*sin(360/n)*0.5
    double r = 1.;
    int n = nbv;
    double computed_area = (r*r)*n*sin(360/n*0.01745329)*0.5;
    INFO("regular polygon area is: " << computed_area);
    CHECK(std::abs(area - computed_area) < 1e-2);
}

TEST_CASE("Test tetra geom", "[geom]") {

    // Construct a regular tetrahedron
    Tetrahedra m;

    m.points.create_points(4);
    m.points[3] = vec3(1/3.0*sqrt(3.0),0,0);
    m.points[2] = vec3(-1/6.0*sqrt(3.0),1/2.0,0);
    m.points[1] = vec3(-1/6.0*sqrt(3.0),-1/2.0,0);
    m.points[0] = vec3(0,0,1/3.0*sqrt(6.0));

    m.create_cells(1);
    m.vert(0, 0) = 0;
    m.vert(0, 1) = 1;
    m.vert(0, 2) = 2;
    m.vert(0, 3) = 3;

    // Check bary
    auto c = m.iter_cells().begin().data;
    auto tet_c = c.geom<Tetrahedron3>();

    vec3 bary{0,0,sqrt(6.0)/12.0};
    INFO("tet bary: " << tet_c.bary_verts());
    CHECK(std::abs((tet_c.bary_verts() - bary).norm2()) < 1e-4);

    // Compute volume of a regular tetrahedron of side length 1 with the formula V = (a^3*sqrt(2))/12
    double v = sqrt(2.0) / 12.0;
    INFO("regular tet volume: " << v);
    INFO("tet volume: " << tet_c.volume());
    CHECK(std::abs(tet_c.volume() - v) < 1e-4);

    // Get first cell facet (match with points at indexes {3,2,1})
    // It's the surface on plane x, y
    m.connect();
    auto f = c.iter_facets().begin().data;
    auto tet_f = f.geom<Triangle3>();
    INFO("facet points: " << tet_f.v[0] << "," << tet_f.v[1] << "," << tet_f.v[2]);

    // Check normal
    INFO("facet normal: " << tet_f.normal());
    CHECK(tet_f.normal().z < 0);
    // Check bary
    INFO("facet bary: " << tet_f.bary_verts());
    CHECK(std::abs((tet_f.bary_verts() - vec3{0,0,0}).norm2()) < 1e-4);
    // Compute area of a facet of a regular tetrahedron (equilateral triangle) of side length 1 with the formula A = (1/4)*sqrt(3)*s²
    double area = 1/4.0*sqrt(3);
    INFO("regular tet facet area: " << area);
    INFO("tet facet area: " << tet_f.unsigned_area());
    CHECK(std::abs(tet_f.unsigned_area() - area) < 1e-4);
}

TEST_CASE("Test hexa geom", "[geom]") {

    // Construct a cube
    Hexahedra m;

    m.points.create_points(8);
    m.points[0] = vec3(-0.5,-0.5,-0.5);
    m.points[1] = vec3(0.5,-0.5,-0.5);
    m.points[2] = vec3(-0.5,0.5,-0.5);
    m.points[3] = vec3(0.5,0.5,-0.5);

    m.points[4] = vec3(-0.5,-0.5,0.5);
    m.points[5] = vec3(0.5,-0.5,0.5);
    m.points[6] = vec3(-0.5,0.5,0.5);
    m.points[7] = vec3(0.5,0.5,0.5);

    m.create_cells(1);
    m.vert(0, 0) = 0;
    m.vert(0, 1) = 1;
    m.vert(0, 2) = 2;
    m.vert(0, 3) = 3;
    m.vert(0, 4) = 4;
    m.vert(0, 5) = 5;
    m.vert(0, 6) = 6;
    m.vert(0, 7) = 7;

    // Check bary
    auto c = m.iter_cells().begin().data;
    auto hex_c = c.geom<Hexahedron3>();

    INFO("hex bary: " << hex_c.bary_verts());
    CHECK(std::abs((hex_c.bary_verts() - vec3{0,0,0}).norm2()) < 1e-4);

    // Compute volume of a hex cube of side length 1 with the formula V = w*h*d
    double v = 1;
    INFO("regular hex volume: " << v);
    INFO("hex volume: " << hex_c.volume());
    CHECK(std::abs(hex_c.volume() - v) < 1e-4);

    // Get first cell facet (match with points at indexes {0,2,4,6})
    // It's the left of the cube
    m.connect();
    auto f = c.iter_facets().begin().data;
    auto hex_f = f.geom<Quad3>();
    INFO("facet points: " << hex_f.v[0] << "," << hex_f.v[1] << "," << hex_f.v[2] << "," << hex_f.v[3]);

    // Check normal
    INFO("facet normal: " << hex_f.normal());
    CHECK(hex_f.normal().x < 0);
    // Check bary
    INFO("facet bary: " << hex_f.bary_verts());
    CHECK(std::abs((hex_f.bary_verts() - vec3{-0.5,0,0}).norm2()) < 1e-4);
    // Compute area of a facet of a cube (square) of side length 1 with the formula A = s²
    double area = 1;
    INFO("regular hex facet area: " << area);
    INFO("hex facet area: " << hex_f.unsigned_area());
    CHECK(std::abs(hex_f.unsigned_area() - area) < 1e-4);
}