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

vec2 rand_v2() {
	return {(rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5};
}

vec3 rand_v3() {
	return {(rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5, (rand()%100 / 100.) - .5};
}

bool are_doubles_equal(double a, double b, double eps) {
	double diff = std::fabs(a - b);
	double max = std::max(std::fabs(a), std::fabs(b));

	// Use absolute tolerance for very small numbers
	if (max < 1e-10) {
		return diff < eps;
	}

	// Otherwise, use relative tolerance
	return diff < max * eps;
}

template<int n>
bool are_vec_equal(vec<n> a, vec<n> b, double eps) {

	vec<n> diff = {std::fabs(a.x - b.x), std::fabs(a.y - b.y)};
	vec<n> max = {std::max(std::fabs(a.x), std::fabs(b.x)), std::max(std::fabs(a.y), std::fabs(b.y))};

	bool equal = true;
	for (int i = 0; i < n; i++) {

		// Check NaN / infinity consistency
		if (std::isnan(a[i]) && std::isnan(b[i]) || std::isinf(a[i]) && std::isinf(b[i]))
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
		auto s = h.geom();

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
	Segment3 s{{sa, sb}};

	INFO("segment: " << s);
	INFO("length: " << s.length());
	
	vec3 p0{2.8, 0, 0.5};
	vec3 p1{-3.8, 0.2, 0.5};
	
	// Check closest point
	INFO("closest point of p0: " << p0 << " is " << s.closest_point(p0));
	INFO("closest point of p1: " << p1 << " is " << s.closest_point(p1));
	CHECK(std::abs((s.closest_point(p0) - sb).norm2()) < 1e-4);
	CHECK(std::abs((s.closest_point(p1) - sa).norm2()) < 1e-4);
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
	Triangle2 t2 = tri_f.project();
	const vec2 &a = t2.v[0];
	const vec2 &b = t2.v[1];
	// a is near from (0,0)
	CHECK(std::abs(a.norm2()) < 1e-4);
	// b is near from the norm of b-a
	vec2 x{sqrt(1), 0};
	CHECK(std::abs((x - b).norm2()) < 1e-4);
	// c is somewhere

	// Check projection xy
	Triangle2 expected_xy_tri{{{tri_f.v[0].x, tri_f.v[0].y}, {tri_f.v[1].x, tri_f.v[1].y}, {tri_f.v[2].x, tri_f.v[2].y}}};
	Triangle2 actual_xy_tri = tri_f.xy();
	INFO("xy projection: " << actual_xy_tri.v[0] << "," << actual_xy_tri.v[1] << "," << actual_xy_tri.v[2]);
	CHECK(std::abs((expected_xy_tri.v[0] - actual_xy_tri.v[0]).norm2()) < 1e-4);
	CHECK(std::abs((expected_xy_tri.v[1] - actual_xy_tri.v[1]).norm2()) < 1e-4);
	CHECK(std::abs((expected_xy_tri.v[2] - actual_xy_tri.v[2]).norm2()) < 1e-4);

	// Check dilate
	auto actual_xy_tri_dilated = actual_xy_tri.dilate(2.);
	CHECK(std::abs(actual_xy_tri.signed_area() - actual_xy_tri_dilated.signed_area() / 4.) < 1e-4);

	// Check unproject xy0
	Triangle3 actual_xy0_tri = actual_xy_tri.xy0();
	INFO("xy0 unproject: " << actual_xy0_tri.v[0] << "," << actual_xy0_tri.v[1] << "," << actual_xy0_tri.v[2]);
	CHECK(std::abs((tri_f.v[0].xy().xy0() - actual_xy0_tri.v[0]).norm2()) < 1e-4);
	CHECK(std::abs((tri_f.v[1].xy().xy0() - actual_xy0_tri.v[1]).norm2()) < 1e-4);
	CHECK(std::abs((tri_f.v[2].xy().xy0() - actual_xy0_tri.v[2]).norm2()) < 1e-4);

	// Check bary coordinates
	vec2 g = actual_xy_tri.bary_verts();
	vec3 expected_bary_coords = vec3{1, 1, 1} / 3.;
	CHECK(std::abs((actual_xy_tri.bary_coords(g) - expected_bary_coords).norm2()) < 1e-4);
	// Check bary on inversed tri
	Triangle3 inversed_t{{{0,0,0},{0.5,0.5,0},{1,0,0}}};
	vec3 g33 = inversed_t.bary_verts();
	CHECK(std::abs((inversed_t.bary_coords(g33) - expected_bary_coords).norm2()) < 1e-4);

	// Check bary coordinates of Triangle3
	vec3 g3 = actual_xy0_tri.bary_verts();
	CHECK(std::abs((actual_xy0_tri.bary_coords(g3) - expected_bary_coords).norm2()) < 1e-4);   

	// Check gradient on some random 2D triangles
	for (int i = 0; i < 10000; i++) {

		Triangle2 r_tri{{rand_v2(), rand_v2(), rand_v2()}};
		vec3 abc = rand_v3();

		vec2 grad = r_tri.grad(abc);

		// Rotate 90° all tri vectors (get orthogonal vectors)
		vec2 o0 = vec2{-r_tri.v[1].y+r_tri.v[2].y, r_tri.v[1].x-r_tri.v[2].x}.normalized();
		vec2 o1 = vec2{-r_tri.v[2].y+r_tri.v[0].y, r_tri.v[2].x-r_tri.v[0].x}.normalized();
		vec2 o2 = vec2{-r_tri.v[0].y+r_tri.v[1].y, r_tri.v[0].x-r_tri.v[1].x}.normalized();
		// Compute lengths of vectors
		double l0 = abc[0] / (o0*(r_tri.v[0] - r_tri.v[1]));
		double l1 = abc[1] / (o1*(r_tri.v[1] - r_tri.v[2]));
		double l2 = abc[2] / (o2*(r_tri.v[2] - r_tri.v[0]));
		// Compute gradient vector
		vec2 exp_grad = (o0 * l0 + o1 * l1 + o2 * l2);

		INFO("check gradient on tri: " << r_tri.v[0] << "," << r_tri.v[1] << "," << r_tri.v[2] << ",");
		INFO("with values :" << abc);
		INFO("grad: " << grad);
		INFO("exp grad: " << exp_grad);
		CHECK(are_vec_equal(grad, exp_grad, 1e-4));
	}

	// TODO Check gradient on some random 3D triangles
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

	// Check quality metrics

	// Results obtained with verdict library formated as following: {quad_coordinate, result}
	VerdictResult<4,2> verdict_jacobian_scale_results[5] {
		{{{0.3,0.6},{1.7,0.5},{1.3,1.5},{0.6,1.2}},0.645942},
		{{{0.9,0.1},{1.2,0.7},{1,1.9},{0.3,1.6}},0.588172},
		{{{0,0.6},{1.2,0.6},{1.1,1.8},{0.7,1.9}},0.880471},
		{{{0.2,0},{1.2,0.3},{1.7,1.5},{0.9,1.2}},0.63186},
		{{{0.2,0.8},{1.9,0.7},{1.3,1.6},{0.1,1.2}},0.798041}
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
		auto custom_quad = custom_m_f.geom<Quad3>().xy();

		// Check consistency between verdict result & ultimaille result
		double scaled_jacobian = custom_quad.scaled_jacobian();
		INFO("jacobian scale [verdict]: " << verdict_result.result);
		INFO("jacobian scale [ultimaille]: " << scaled_jacobian);
		CHECK(std::abs(scaled_jacobian - verdict_result.result) < 1e-5);
	}

	// Check signed area of Quad2
	// Create a rect
	Quad2 q2{{{0,0}, {2,0}, {2,1}, {0,1}}};
	double actual_q2_signed_area = q2.signed_area();
	// Using the rect formula of area: a=l*h
	double expected_signed_area = (q2.v[1] - q2.v[0]).norm() * (q2.v[3] - q2.v[0]).norm();
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

	// Check bary coordinates
	vec3 tet_c_g = tet_c.bary_verts();
	vec4 expected_tet_bary_coords = vec4{1, 1, 1, 1} / 4.;
	CHECK(std::abs((tet_c.bary_coords(tet_c_g) - expected_tet_bary_coords).norm2()) < 1e-4);

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
	Triangle3 tri{{{0,0,0},{1,0,0},{0.5,0.5,0}}};
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

	// Check quality metrics

	// Results obtained with verdict library formated as following: {quad_coordinate, result}
	VerdictResult<8,3> verdict_jacobian_scale_results[10] {
		{{{-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.5,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},1},
		{{{-0.4,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.6,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.977069},
		{{{-0.3,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.7,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.916841},
		{{{-0.2,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.8,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.832049},
		{{{-0.1,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.9,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.734169},
		{{{0,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.632456},
		{{{0.1,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.1,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.53357},
		{{{0.2,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.2,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.441714},
		{{{0.3,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.3,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.359086},
		{{{0.4,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,1.4,-0.5},{0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}},0.286442}
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
		auto custom_quad = custom_m_f.geom<Hexahedron3>();

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
	auto pyr_c = c.geom<Pyramid3>();

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