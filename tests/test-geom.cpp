#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>
#include <chrono>
#include <cmath>
#include <algorithm>

using namespace UM;

// This structure will contains some results of verdict quality library to compare
// https://github.com/sandialabs/verdict
template<int n, int m>
struct VerdictResult {
    double coordinates[n][m];
    double result;
};

// Some useful functions

double rand_scalar() {
    return rand()%100 / 100.;
}

double rand_scalar(float from, float to) {
    return from + rand()%100 / 100. * (to - from);
}

vec2 rand_v2() {
    return {(rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5};
}

vec3 rand_v3() {
    return {(rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5};
}

vec4 rand_v4() {
    return {(rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5};
}

bool are_doubles_equal(double a, double b, double eps) {
    double diff = std::fabs(a - b);
    double max = std::max(std::fabs(a), std::fabs(b));

    // Check NaN / infinity consistency
    if ((std::isnan(a) && std::isnan(b)) || (std::isinf(a) && std::isinf(b)))
        return true;

    // Use absolute tolerance for very small numbers
    if (max < 1e-10) {
        return diff < eps;
    }

    // Otherwise, use relative tolerance
    return diff < max * eps;
}

template<int n>
bool are_vec_equal(vec<n> a, vec<n> b, double eps) {

    vec<n> diff;
    vec<n> max;
    for (int i = 0; i < n; i++) {
        diff[i] = std::fabs(a[i] - b[i]);
        max[i] = std::max(std::fabs(a[i]), std::fabs(b[i]));
    }

    // vec<n> diff = {std::fabs(a.x - b.x), std::fabs(a.y - b.y)};
    // vec<n> max = {std::max(std::fabs(a.x), std::fabs(b.x)), std::max(std::fabs(a.y), std::fabs(b.y))};

    bool equal = true;
    for (int i = 0; i < n; i++) {

        // Check NaN / infinity consistency
        if ((std::isnan(a[i]) && std::isnan(b[i])) || (std::isinf(a[i]) && std::isinf(b[i])))
            continue;

        // Use absolute tolerance for very small numbers
        if (max[i] < 1e-10)
            equal &= diff[i] < eps;
        else 
            // Otherwise, use relative tolerance
            equal &= diff[i] < max[i] * eps;
    }

    return equal;
}

TEST_CASE("aggregate", "[geom]") {

    CHECK(std::is_aggregate<Segment2>());
    CHECK(std::is_aggregate<Segment3>());
    // CHECK(std::is_aggregate<Triangle2>());
    // CHECK(std::is_aggregate<Triangle3>());
    // CHECK(std::is_aggregate<Quad2>());
    // CHECK(std::is_aggregate<Quad3>());
    // CHECK(std::is_aggregate<Poly3>());
    // CHECK(std::is_aggregate<Tetrahedron>());
    // CHECK(std::is_aggregate<Hexahedron>());
    // CHECK(std::is_aggregate<Pyramid>());

    CHECK(std::is_aggregate<vec2>());
    CHECK(std::is_aggregate<vec3>());
}

// TEST_CASE("vec benchmark", "[geom]") {

//     auto start = std::chrono::high_resolution_clock::now();
// 	for (int i = 0; i < 1000000000; i++) {
// 		vec2 v{rand()%100, rand()%100};
// 	}
//     // Capture the end time
//     auto end = std::chrono::high_resolution_clock::now();

//     // Compute the duration
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

//     // Output the duration
//     INFO("DURATION: " << duration.count() << " milliseconds");
// 	CHECK(false);
// }

TEST_CASE("Test segment geom", "[geom]") {

    // Create mesh with one equilateral triangle cell
    Triangles m;
    m.points.create_points(3);
    m.points[0] = {0,0,0};
    m.points[1] = vec3(1,0,0);
    m.points[2] = vec3(0.5,0.8660,0);
    m.create_facets(1);
    m.vert(0, 0) = 0;
    m.vert(0, 1) = 1;
    m.vert(0, 2) = 2;
    m.connect();

    // Get first half-edge and extract segment
    for (int hi = 0; hi < 3; hi++) {
        auto h = Surface::Halfedge(m, hi);
        Segment3 s = h;

        // Check that segment references points correctly
        CHECK(std::abs((s[0] - m.points[hi % 3]).norm2()) < 1e-4);
        CHECK(std::abs((s[1] - m.points[(hi + 1) % 3]).norm2()) < 1e-4);

        // Check length
        INFO("halfedge " << hi << " - length of segment: " << s.length());
        CHECK(std::abs(s.length() - 1.) < 1e-4);
    }

    // Create a segment
    vec3 sa{0,0,0};
    vec3 sb{1,0,0};
    Segment3 s(sa, sb);

    INFO("segment: " << s);
    INFO("length: " << s.length());

    vec3 p0{2.8, 0, 0.5};
    vec3 p1{-3.8, 0.2, 0.5};

    // Check closest point
    INFO("closest point of p0: " << p0 << " is " << s.nearest_point(p0));
    INFO("closest point of p1: " << p1 << " is " << s.nearest_point(p1));
    CHECK(std::abs((s.nearest_point(p0) - sb).norm2()) < 1e-4);
    CHECK(std::abs((s.nearest_point(p1) - sa).norm2()) < 1e-4);
}

TEST_CASE("Test triangle3 nearest point", "[geom]") {
    Triangle3 t = {{-1, 5, 0}, {2, 2, -3}, {5, 5, 0}};
    std::vector<vec3> pts = {{1, 1, 1}, {-1, -3, -4}, {2, 4, -1}, {-2.732051, 6.732051, 1.732051}, {3, 7, -4}};
    std::vector<vec3> prj = {{1, 3.5, -1.5}, {2, 2, -3}, {2, 4, -1}, {-1, 5, 0}, {3, 4, -0.9999995}};
    for (auto [p,q] : zip(pts, prj)) {
        vec3 nearest = t.nearest_point(p);
        CHECK( (nearest-q).norm() < 1e-4 );
    }
}


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
        CHECK(Triangle3(f).normal().z > 0);
    }

    // Check bary
    auto f = m.iter_facets().begin().f;
    Triangle3 tri_f = f;
    CHECK(std::abs((tri_f.bary_verts() - vec3{0.5,0.5/3.,0}).norm()) < 1e-4);

    // Check area comparing with formula A = bh / 2
    double tri_area = tri_f.unsigned_area();
    double test_area = 1. /* b */ * 0.5 /* h */ * 0.5 /* divide by two */;
    CHECK(std::abs(tri_area - test_area) < 1e-4);

    // Check projection
    Triangle2 t2 = tri_f.project();
    const vec2 &a = t2[0];
    const vec2 &b = t2[1];
    // a is near from (0,0)
    CHECK(std::abs(a.norm2()) < 1e-4);
    // b is near from the norm of b-a
    vec2 x{sqrt(1), 0};
    CHECK(std::abs((x - b).norm2()) < 1e-4);
    // c is somewhere

    // Check projection xy
    Triangle2 expected_xy_tri({tri_f[0].x, tri_f[0].y}, {tri_f[1].x, tri_f[1].y}, {tri_f[2].x, tri_f[2].y});
    Triangle2 actual_xy_tri = tri_f.xy();
    INFO("xy projection: " << actual_xy_tri[0] << "," << actual_xy_tri[1] << "," << actual_xy_tri[2]);
    CHECK(std::abs((expected_xy_tri[0] - actual_xy_tri[0]).norm2()) < 1e-4);
    CHECK(std::abs((expected_xy_tri[1] - actual_xy_tri[1]).norm2()) < 1e-4);
    CHECK(std::abs((expected_xy_tri[2] - actual_xy_tri[2]).norm2()) < 1e-4);

    // Check dilate
    auto actual_xy_tri_dilated = actual_xy_tri.dilate(2.);
    CHECK(std::abs(actual_xy_tri.signed_area() - actual_xy_tri_dilated.signed_area() / 4.) < 1e-4);

    // Check unproject xy0
    Triangle3 actual_xy0_tri = actual_xy_tri.xy0();
    INFO("xy0 unproject: " << actual_xy0_tri[0] << "," << actual_xy0_tri[1] << "," << actual_xy0_tri[2]);
    CHECK(std::abs((tri_f[0].xy().xy0() - actual_xy0_tri[0]).norm2()) < 1e-4);
    CHECK(std::abs((tri_f[1].xy().xy0() - actual_xy0_tri[1]).norm2()) < 1e-4);
    CHECK(std::abs((tri_f[2].xy().xy0() - actual_xy0_tri[2]).norm2()) < 1e-4);

    // Check bary coordinates
    vec2 g = actual_xy_tri.bary_verts();
    vec3 expected_bary_coords = vec3{1, 1, 1} / 3.;
    CHECK(std::abs((actual_xy_tri.bary_coords(g) - expected_bary_coords).norm2()) < 1e-4);
    // Check bary on inversed tri
    Triangle3 inversed_t({0,0,0},{0.5,0.5,0},{1,0,0});
    vec3 g33 = inversed_t.bary_verts();
    CHECK(std::abs((inversed_t.bary_coords(g33) - expected_bary_coords).norm2()) < 1e-4);

    // Check bary coordinates of Triangle3
    vec3 g3 = actual_xy0_tri.bary_verts();
    CHECK(std::abs((actual_xy0_tri.bary_coords(g3) - expected_bary_coords).norm2()) < 1e-4);   

    // Check gradient on some random 2D triangles
    for (int i = 0; i < 10000; i++) {

        Triangle2 r_tri(rand_v2(), rand_v2(), rand_v2());
        vec3 abc = rand_v3();

        vec2 grad = r_tri.grad(abc);

        // Rotate 90° all tri vectors (get orthogonal vectors)
        vec2 o0 = vec2{-r_tri[1].y+r_tri[2].y, r_tri[1].x-r_tri[2].x}.normalized();
        vec2 o1 = vec2{-r_tri[2].y+r_tri[0].y, r_tri[2].x-r_tri[0].x}.normalized();
        vec2 o2 = vec2{-r_tri[0].y+r_tri[1].y, r_tri[0].x-r_tri[1].x}.normalized();
        // Compute lengths of vectors
        double l0 = abc[0] / (o0*(r_tri[0] - r_tri[1]));
        double l1 = abc[1] / (o1*(r_tri[1] - r_tri[2]));
        double l2 = abc[2] / (o2*(r_tri[2] - r_tri[0]));
        // Compute gradient vector
        vec2 exp_grad = (o0 * l0 + o1 * l1 + o2 * l2);

        INFO("check gradient on tri: " << r_tri[0] << "," << r_tri[1] << "," << r_tri[2] << ",");
        INFO("with values :" << abc);
        INFO("grad: " << grad);
        INFO("exp grad: " << exp_grad);
        CHECK(are_vec_equal(grad, exp_grad, 1e-4));

        // for (int v = 0; v < 3; v++) {
        // 	INFO("grad:" << grad);
        // 	INFO("a:" << grad * (r_tri[(v + 1) % 3] - r_tri[v]));
        // 	INFO("b:" << abc[(v + 1) % 3] - abc[v]);
        // 	CHECK(are_doubles_equal(grad * (r_tri[(v + 1) % 3] - r_tri[v]), abc[(v + 1) % 3] - abc[v], 1e-4));
        // }
    }

    // // Check gradient on some random 3D triangles
    // for (int i = 0; i < 10000; i++) {

    // 	Triangle3 r_tri(rand_v3(), rand_v3(), rand_v3());
    // 	vec3 abc = rand_v3();

    // 	vec3 grad = r_tri.grad(abc);

    // 	// Get orthogonal vectors
    // 	vec3 o0 = cross(cross(r_tri[1] - r_tri[0], r_tri[2] - r_tri[0]), r_tri[2] - r_tri[1])/*.normalized()*/;
    // 	vec3 o1 = cross(cross(r_tri[1] - r_tri[0], r_tri[2] - r_tri[1]), r_tri[2] - r_tri[0]).normalized();
    // 	vec3 o2 = cross(cross(r_tri[0] - r_tri[2], r_tri[1] - r_tri[2]), r_tri[1] - r_tri[0]).normalized();

    // 	// Compute lengths of vectors
    // 	double l0 = abc[0] / (o0*(r_tri[0] - r_tri[1]));
    // 	double l1 = abc[1] / (o1*(r_tri[1] - r_tri[2]));
    // 	double l2 = abc[2] / (o2*(r_tri[2] - r_tri[0]));
    // 	// Compute gradient vector
    // 	vec3 exp_grad = (o0 * l0 + o1 * l1 + o2 * l2);

    // 	INFO("check gradient on tri: " << r_tri[0] << "," << r_tri[1] << "," << r_tri[2] << ",");
    // 	INFO("with values :" << abc);
    // 	INFO("grad: " << grad);
    // 	INFO("exp grad: " << exp_grad);
    // 	CHECK(are_vec_equal(grad, exp_grad, 1e-4));
    // 	// CHECK(false);
    // }

    // TODO Check gradient on some random tetrahedrons

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
    Quad3 quad_f = f;

    INFO("quad: " << quad_f.normal());

    // Check normal
    for (auto f : m.iter_facets()) {
        auto n = Quad3(f).normal().z;
        CHECK(n > 0);
    }

    // Check area
    double quad_area = quad_f.unsigned_area();
    // Considering the quad is a square of length 1 and height 1
    // We have an area equal to 1*1=1
    bool is_area_correct = quad_area >= 0.9999 && quad_area <= 1.0001;
    CHECK(is_area_correct);

    // Check quality metrics

    // Results obtained with verdict library formated as following: {quad_coordinate, result}
    VerdictResult<4,2> verdict_jacobian_scale_results[6] {
        {{{0.3,0.6},{1.7,0.5},{1.3,1.5},{0.6,1.2}},0.645942},
        {{{0.9,0.1},{1.2,0.7},{1,1.9},{0.3,1.6}},0.588172},
        {{{0,0.6},{1.2,0.6},{1.1,1.8},{0.7,1.9}},0.880471},
        {{{0.2,0},{1.2,0.3},{1.7,1.5},{0.9,1.2}},0.63186},
        {{{0.2,0.8},{1.9,0.7},{1.3,1.6},{0.1,1.2}},0.798041},
        {{{0,0},{1e-31,0},{1,1},{0,1}},0}
    };

    // Construct quads from verdict results
    for (auto verdict_result : verdict_jacobian_scale_results) {

        // Construct a quad from given coordinates
        Quads custom_m;
        custom_m.points.create_points(6);
        for (int i = 0; i < 4; i++) {
            custom_m.points[i] = vec3(verdict_result.coordinates[i][0], verdict_result.coordinates[i][1], 0);
        }

        custom_m.create_facets(1);
        custom_m.vert(0, 0) = 0;
        custom_m.vert(0, 1) = 1;
        custom_m.vert(0, 2) = 2;
        custom_m.vert(0, 3) = 3;

        // Get geometry of first face
        auto custom_m_f = custom_m.iter_facets().begin().f;
        auto custom_quad = Quad3(custom_m_f).xy();

        // Check consistency between verdict result & ultimaille result
        double scaled_jacobian = custom_quad.scaled_jacobian();
        INFO("jacobian scale [verdict]: " << verdict_result.result);
        INFO("jacobian scale [ultimaille]: " << scaled_jacobian);
        CHECK(std::abs(scaled_jacobian - verdict_result.result) < 1e-5);
    }

    // Check signed area of Quad2
    // Create a rect
    Quad2 q2({0,0}, {2,0}, {2,1}, {0,1});
    double actual_q2_signed_area = q2.signed_area();
    // Using the rect formula of area: a=l*h
    double expected_signed_area = (q2[1] - q2[0]).norm() * (q2[3] - q2[0]).norm();
    INFO("signed area of q2 is: " << actual_q2_signed_area);
    CHECK(std::abs(actual_q2_signed_area - expected_signed_area) < 1e-4);

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
    Poly3 poly_f = f;

    // // Check normal
    for (auto f : m.iter_facets()) {
        CHECK(Poly3(f).normal().z > 0);
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
    Tetrahedron tet_c = c;

    vec3 bary{0,0,sqrt(6.0)/12.0};
    INFO("tet bary: " << tet_c.bary_verts());
    CHECK(std::abs((tet_c.bary_verts() - bary).norm2()) < 1e-4);

    // Compute volume of a regular tetrahedron of side length 1 with the formula V = (a^3*sqrt(2))/12
    double v = sqrt(2.0) / 12.0;
    INFO("regular tet volume: " << v);
    INFO("tet volume: " << tet_c.volume());
    CHECK(std::abs(tet_c.volume() - v) < 1e-4);

    // // Check bary coordinates
    // {
    // 	vec3 tet_c_g = tet_c.bary_verts();
    // 	vec4 expected_tet_bary_coords = vec4{1, 1, 1, 1} / 4.;
    // 	CHECK(std::abs((tet_c.bary_coords(tet_c_g) - expected_tet_bary_coords).norm2()) < 1e-4);
    // }

    {

        for (int i = 0; i < 10000; i++) {
            // Compute random vec4 that sum up to 1 !
            float a = 1.;
            float u = rand_scalar(0., 1.);
            a -= u;
            float v = rand_scalar(0., a);
            a -= v;
            float w = rand_scalar(0., a);
            float w2 =  1. - (u + v + w);

            // vec4 expected_bary_coords{0.2, 0.4, 0.15, 0.25};
            vec4 expected_bary_coords{u, v, w, w2};
            // Barycentric coordinates must sum up to 1 (here this checking shouldn't be necessary because of computation of u,v,w,w' above)
            REQUIRE(std::abs(expected_bary_coords[0] + expected_bary_coords[1] + expected_bary_coords[2] + expected_bary_coords[3] - 1) < 1e-4);

            // Compute point from barycentric coordinates
            // Doesn't want a flat tet with volume = 0
            vec3 A = vec3{0,0,0} + rand_v3() * .5;
            vec3 B = vec3{1,0,0} + rand_v3() * .5;
            vec3 C = vec3{0,1,0} + rand_v3() * .5;
            vec3 D = vec3{1,1,0} + rand_v3() * .5;
            // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates.html
            // u, v, w, w' the barycentric coordinates
            // u = 0.2, v = 0.4, w = 0.15, w' = 0.25
            // P = A + u * AB + v * AC + w * AD 
            // so w'= (1 - u - v - w) then P = (1 - u - v - w) * A + u * B + v * C + w * D
            // match with P = w' * A + u * B + v * C + w * D instead of P = u * A + v * B + w * C + w' * D
            // so I shift bary coordinates to get the right order in coordinates components...
            vec3 P =  A + expected_bary_coords[1] * (B - A) + expected_bary_coords[2] * (C - A) + expected_bary_coords[3] * (D - A);

            // Compute barycentric coordinates from point
            Tetrahedron tet_c2(A, B, C, D);
            vec4 actual_bary_coords = tet_c2.bary_coords(P);

            CHECK(std::abs((expected_bary_coords - actual_bary_coords).norm2()) < 1e-4);
        }
    }

    // Get first cell facet (match with points at indexes {3,2,1})
    // It's the surface on plane x, y
    m.connect();
    auto f = c.iter_facets().begin().data;
    Triangle3 tet_f = f;
    INFO("facet points: " << tet_f[0] << "," << tet_f[1] << "," << tet_f[2]);

    // Check normal
    INFO("facet normal: " << tet_f.normal());
    CHECK(tet_f.normal().z < 0);
    // Check bary
    INFO("facet bary: " << tet_f.bary_verts());
    CHECK(std::abs((tet_f.bary_verts() - vec3{0,0,0}).norm2()) < 1e-4);
    // Check bary coordinates
    vec3 g{0,0,0};
    INFO("facet bary coords: " << tet_f.bary_coords(g));
    // As it's a regular tet, with regular tri as facets we expect the following bary coords
    vec3 expected_bary_coords = vec3{1, 1, 1} / 3.;
    CHECK(std::abs((tet_f.bary_coords(g) - expected_bary_coords).norm2()) < 1e-4);
    // Compute area of a facet of a regular tetrahedron (equilateral triangle) of side length 1 with the formula A = (1/4)*sqrt(3)*s²
    double area = 1/4.0*sqrt(3);
    INFO("regular tet facet area: " << area);
    INFO("tet facet area: " << tet_f.unsigned_area());
    CHECK(std::abs(tet_f.unsigned_area() - area) < 1e-4);
    // Check dilate
    auto tet_f_dilated = tet_f.dilate(2.);
    CHECK(std::abs(tet_f.unsigned_area() - tet_f_dilated.unsigned_area() / 4.) < 1e-4);
    // Check tangent basis
    mat3x3 t = tet_f.tangent_basis({1,0,0});
    INFO("tangent basis: " << t);
    CHECK(std::abs((t[1] - vec3{0,-1,0}).norm2()) < 1e-4);
    // Check tagent basis
    Triangle3 tri({0,0,0},{1,0,0},{0.5,0.5,0});
    t = tri.tangent_basis();
    INFO("tangent basis: " << t);
    CHECK(std::abs((t[1] - vec3{0,1,0}).norm2()) < 1e-4);
    CHECK(std::abs((t[2] - vec3{0,0,1}).norm2()) < 1e-4);
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
    Hexahedron hex_c = c;

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
    Quad3 hex_f = f;
    INFO("facet points: " << hex_f[0] << "," << hex_f[1] << "," << hex_f[2] << "," << hex_f[3]);

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

    // Check quality metrics

    // Results obtained with verdict library formated as following: {quad_coordinate, result}
    VerdictResult<8,3> verdict_jacobian_scale_results[11] {
        {{{-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.5,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},1},
        {{{-0.4,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.6,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.977069},
        {{{-0.3,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.7,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.916841},
        {{{-0.2,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.8,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.832049},
        {{{-0.1,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.9,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.734169},
        {{{0,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.632456},
        {{{0.1,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.1,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.53357},
        {{{0.2,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.2,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.441714},
        {{{0.3,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.3,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.359086},
        {{{0.4,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.4,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.286442},
        {{{0,0,0},{0,0,0},{-0.5,0.5,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},1e+30}
    };

    // Construct quads from verdict results
    for (auto verdict_result : verdict_jacobian_scale_results) {

        // Construct a quad from given coordinates
        Hexahedra custom_m;
        custom_m.points.create_points(8);
        for (int i = 0; i < 8; i++) {
            custom_m.points[i] = vec3(verdict_result.coordinates[i][0], verdict_result.coordinates[i][1], verdict_result.coordinates[i][2]);
        }

        custom_m.create_cells(1);
        custom_m.vert(0, 0) = 0;
        custom_m.vert(0, 1) = 1;
        custom_m.vert(0, 2) = 2;
        custom_m.vert(0, 3) = 3;
        custom_m.vert(0, 4) = 4;
        custom_m.vert(0, 5) = 5;
        custom_m.vert(0, 6) = 6;
        custom_m.vert(0, 7) = 7;

        // Get geometry of first face
        auto custom_m_f = custom_m.iter_cells().begin().data;
        Hexahedron custom_quad = custom_m_f;

        // Check consistency between verdict result & ultimaille result
        double scaled_jacobian = custom_quad.scaled_jacobian();
        INFO("jacobian scale [verdict]: " << verdict_result.result);
        INFO("jacobian scale [ultimaille]: " << scaled_jacobian);
        CHECK(std::abs(scaled_jacobian - verdict_result.result) < 1e-5);
    }

}

TEST_CASE("Test pyramid geom", "[geom]") {

    // Construct an equilateral square pyramid
    // https://en.wikipedia.org/wiki/Square_pyramid
    Pyramids m;

    m.points.create_points(5);
    double h = sqrt(0.5);
    // Base
    m.points[0] = vec3(-0.5,-0.5,0);
    m.points[1] = vec3(0.5,-0.5,0);
    m.points[2] = vec3(0.5,0.5,0);
    m.points[3] = vec3(-0.5,0.5,0);
    // Apex
    m.points[4] = vec3(0,0,h);

    m.create_cells(1);
    m.vert(0, 0) = 0;
    m.vert(0, 1) = 1;
    m.vert(0, 2) = 2;
    m.vert(0, 3) = 3;
    m.vert(0, 4) = 4;

    // Check bary
    auto c = m.iter_cells().begin().data;
    Pyramid pyr_c = c;

    INFO("pyramid bary: " << pyr_c.bary_verts());
    // All distances from apex should be equal to 1
    INFO("dist a:" << (m.points[0] - m.points[4]).norm());
    INFO("dist b:" << (m.points[1] - m.points[4]).norm());
    INFO("dist c:" << (m.points[2] - m.points[4]).norm());
    INFO("dist d:" << (m.points[3] - m.points[4]).norm());
    CHECK(std::abs((pyr_c.bary_verts() - vec3{0,0,h/5}).norm2()) < 1e-4);

    // Compute volume of a square pyramid of side length 1 with the formula V = 1/3*l²*h
    double v = 1/3.0*h;
    INFO("square pyramid volume: " << v);
    INFO("pyramid volume: " << pyr_c.volume());
    CHECK(std::abs(pyr_c.volume() - v) < 1e-4);

    // Get pyramid base geometry
    // It should be oriented up (Z-up)
    auto base_geom = pyr_c.base();
    // Check base normal
    INFO("base normal: " << base_geom.normal());
    CHECK(base_geom.normal().z > 0);
    // Check base bary
    INFO("base bary: " << base_geom.bary_verts());
    CHECK(std::abs((base_geom.bary_verts() - vec3{0,0,0}).norm2()) < 1e-4);
    // Check base area
    double area = 1;
    INFO("square base area: " << area);
    INFO("base area: " << base_geom.unsigned_area());
    CHECK(std::abs(base_geom.unsigned_area() - area) < 1e-4);
}

TEST_CASE("Test wedge geom", "[geom]") {

    // Random wedge sample
    for (int i = 0; i < 10000; i++) {

        // Check volume

        // Random parameters
        double a = rand_scalar(), b = rand_scalar(), h = rand_scalar();
        // Random start position
        vec3 p_offset = rand_v3();

        // Create random wedge with a rectangular base
        Wedge w(
            p_offset, 
            p_offset + vec3{0,b,0}, 
            p_offset + vec3{a*.25/*rand_scalar()*/,b*.5,h}, 
            p_offset + vec3{a,0,0}, 
            p_offset + vec3{a,b,0}, 
            p_offset + vec3{a-a*.25/*rand_scalar()*/,b*.5,h}
        );
        
        // Compute volume of a wedge with a rectangular base with the formula: V = bh((a/3)+(c/6))
        // ref: https://en.wikipedia.org/wiki/Wedge_(geometry)
        double c = (w[2] - w[5]).norm();
        double exp_v = b*h*((a/3.)+(c/6.));

        double act_v = w.volume();
        
        INFO("expected volume: " << exp_v);
        INFO("actual volume: " << act_v);
        CHECK(std::abs(act_v - exp_v) < 1e-4);

        // Check bary
        vec3 act_bary = w.bary_verts();
        vec3 exp_bary = (w[0] + w[1] + w[2] + w[3] + w[4] + w[5]) / 6.;
        INFO("expected bary: " << exp_bary);
        INFO("actual bary: " << act_bary);
        CHECK(std::abs((act_bary - exp_bary).norm2()) < 1e-4);
    }
}

// This test aims to check that extracting an abstract polygon geometry from a facet 
// give equivalent computations than computation from facet real shape
TEST_CASE("Test polygon facet extraction geom", "[geom]") {
    
    Wedges w;
    w.points.create_points(6);
    w.create_cells(1);
    w.points[0] = {0,0,0};
    w.points[1] = {1,0,0};
    w.points[2] = {0.5,0,1};
    w.points[3] = {0,1,0};
    w.points[4] = {1,1,0};
    w.points[5] = {0.5,1,1};

    for (int i = 0; i < 6; i++)
        w.vert(0, i) = i;

    // Should be a triangle
    CHECK(Poly3(Volume::Facet(w, 0)).v.size() == 3);
    // Should be a quad
    CHECK(Poly3(Volume::Facet(w, 2)).v.size() == 4);


    Pyramids p;
    p.points.create_points(5);
    p.create_cells(5);
    p.points[0] = {0,0,0};
    p.points[1] = {1,0,0};
    p.points[2] = {1,1,0};
    p.points[3] = {0,1,0};
    p.points[4] = {0.5,0.5,1};

    for (int i = 0; i < 5; i++)
        w.vert(0, i) = i;

    // Should be a quad
    CHECK(Poly3(Volume::Facet(p, 0)).v.size() == 4);
    // Should be triangles
    for (int i = 1; i < 5; i++)
        CHECK(Poly3(Volume::Facet(p, i)).v.size() == 3);


}

TEST_CASE("Test segment 2", "[geom][segment]") {
    for (int i = 0; i < 100; i++) {
        vec2 a{0,0};
        vec2 b{rand_scalar(), 0};

        Segment2 s{a, b};
        REQUIRE(std::abs((s.a - a).norm2()) < 1e-4);

        CHECK(std::abs(s.length() - (b.x - a.x)) < 1e-4);
        CHECK(std::abs(s.length() * s.length() - s.length2()) < 1e-4);

        CHECK(std::abs((s.vector() - (b - a)).norm2()) < 1e-4);

        Segment3 s3{{a.x, a.y, 0}, {b.x, b.y, 0}};
        CHECK(std::abs((s.xy0().vector() - s3.vector()).norm2()) < 1e-4);

        // Case: point between segment bounds
        {
            // A point between a and b on x coords and random y coords
            // Should be at distance of y from segment
            double d = rand_scalar();
            vec2 p{rand_scalar(0, b.x), d};
            CHECK(std::abs(s.distance(p) - d) < 1e-4);
        }

        // Case: point out of segment bounds
        {
            double d = rand_scalar();
            vec2 p{rand_scalar(b.x, b.x + d), b.y};
            CHECK(std::abs(s.distance(p) - (p.x - b.x)) < 1e-4);
        }
        {
            double d = rand_scalar();
            vec2 p{rand_scalar(-d, 0), b.y};
            CHECK(std::abs(s.distance(p) - (-p.x)) < 1e-4);
        }

        // Case: degenerated segment
        {
            // Distance of degenerated segment to a point 
            // is distance from point to point
            vec2 p = rand_v2();
            vec2 a = rand_v2();
            vec2 b = a + vec2{0, 1e-16};
            Segment2 degenerated_s{a, b};
            vec2 v = p - a;
            double dist = sqrt(v * v);

            double actual_dist = degenerated_s.distance(p);

            INFO("expected dist: " << dist << ", from: " << degenerated_s << ", to: " << p);
            INFO("actual dist: " << actual_dist);
            CHECK(std::abs(actual_dist - dist) < 1e-4);
        }
    }
}

TEST_CASE("Test segment 3", "[geom][segment]") {
    for (int i = 0; i < 100; i++) {
        vec3 a{0,0,0};
        vec3 b{rand_scalar(), 0, 0};

        Segment3 s{a, b};
        REQUIRE(std::abs((s.a - a).norm2()) < 1e-4);

        CHECK(std::abs(s.length() - (b.x - a.x)) < 1e-4);
        CHECK(std::abs(s.length() * s.length() - s.length2()) < 1e-4);

        CHECK(std::abs((s.vector() - (b - a)).norm2()) < 1e-4);

        Segment2 s2{{a.x, a.y}, {b.x, b.y}};
        CHECK(std::abs((s.xy().vector() - s2.vector()).norm2()) < 1e-4);
    }

    for (int i = 0; i < 1000; i++) {
        // Compute random u, v that sum up to 1 !
        float u = rand_scalar(-2., 2.);
        float v =  1. - u;
        vec2 expected_bary_coords{u, v};
        
        vec3 A = vec3{0,0,0} + rand_v3() * .5;
        vec3 B = vec3{1,0,0} + rand_v3() * .5;
        // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates.html
        // u, v the barycentric coordinates
        // P = A + u * AB
        // so v = (1 - u) then P = (1 - u) * A + u * B
        Segment3 s{A, B};
        vec3 P =  A + v * (B - A);

        double actual_u = s.bary_coords(P);
        INFO("u:" << actual_u << ", " << u << "," << v);

        CHECK(std::abs((actual_u - u) < 1e-4));

        vec3 actual_closest = s.closest_point(P);
        vec3 expected_closest = s.a;
        if (u < 0.)
            expected_closest = s.b;
        else if (u > 1)
            expected_closest = s.a;
        else 
            expected_closest = s.b + u * (s.a - s.b);
            
        CHECK(std::abs((expected_closest - actual_closest).norm2()) < 1e-4);
    }
}
