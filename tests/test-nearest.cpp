#include <catch2/catch_test_macros.hpp>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("projection to a 3d triangulated surface", "[nearest]") {
    Triangles m; // a simple tetrahedron
    *m.points.data = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
    m.facets = {1,3,2, 3,1,0, 2,0,1, 0,2,3};

    NearestPointOnMesh npm(m);

    PointOnMesh p1 = npm.query({0.8311, 0.0224465, 1.00753});
    CHECK( (p1.f == 2) || (p1.f == 3) ); // nearest point is on an edge between facets 2 and 3
    CHECK( (vec3{0.411785, 0, 0.588215} - static_cast<vec3>(p1)).norm() < 1e-5 );

    auto p2 = npm.query({0.802483, 1.14223, 1.26074});
    CHECK( p2.f == 2 ); // nearest point is inside the facet 2
    CHECK( (vec3{0.067334, 0.407079, 0.525587} - static_cast<vec3>(p2)).norm() < 1e-5 );

    auto p3 = npm.query({1.2124, -0.118379, 0.194748});
    CHECK( (p3.f == 1) || (p3.f == 2) || (p3.f == 3) ); // nearest is on the vertex between facets 1, 2 and 3
    CHECK( (vec3{1, 0, 0} - static_cast<vec3>(p1)).norm() < 1e-5 );
}

