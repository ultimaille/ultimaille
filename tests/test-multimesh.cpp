#include <catch2/catch_test_macros.hpp>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("Shared pointset editing", "[MultiMesh]") {
    Triangles m1;
    *m1.points.data = {{1, 0, 0}, { 0.309017, 0.951057, 0}, {-0.809017, 0.587785, 0}};
    m1.facets = {0,1,2};

    PointAttribute<int> a1(m1.points);
    for (int i=0; i<3; i++)
        a1[i] = i;

    REQUIRE( a1.ptr->data.size() == 3 );

    Quads m2;
    m2.points = m1.points;

    PointAttribute<int> a2(m2.points);
    m2.points.push_back({-0.809017,-0.587785,0});
    m2.points.push_back({ 0.309017,-0.951057,0});
    m2.facets = {0,2,3,4};

    REQUIRE( a1.ptr->data.size() == 5 );
    REQUIRE( a2.ptr->data.size() == 5 );

    for (int i=3; i<5; i++)
        a2[i] = i;

    write_by_extension("m1.geogram", m1, {{"a1", a1}, {"a2", a2}});
    write_by_extension("m2.geogram", m2, {{"a1", a1}, {"a2", a2}});
}
