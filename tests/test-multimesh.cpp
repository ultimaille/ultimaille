#include <catch2/catch_test_macros.hpp>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("Shared pointset editing", "[MultiMesh]") {
    Triangles m1;
    *m1.points.data = {{-1,-1,-1}, {1, 0, 0}, { 0.309017, 0.951057, 0}, {-0.809017, 0.587785, 0}};
    m1.facets = {1,2,3};

    write_by_extension("m1-.geogram", m1);

    PointAttribute<int> a1(m1.points);
    for (int i=0; i<4; i++)
        a1[i] = i;

    REQUIRE( a1.ptr->data.size() == 4 );

    Quads m2;
    m2.points = m1.points;

    PointAttribute<int> a2(m2.points);
    m2.points.push_back({-0.809017,-0.587785,0});
    m2.points.push_back({ 0.309017,-0.951057,0});
    m2.facets = {1,3,4,5};

    REQUIRE( a1.ptr->data.size() == 6 );
    REQUIRE( a2.ptr->data.size() == 6 );

    for (int i=4; i<6; i++)
        a2[i] = i;

    std::vector<int> old2new;
    m2.points.delete_points([&](int id) { return id==0; }, old2new);

    for (auto v : range(3))
        m1.vert(0, v) = old2new[m1.vert(0, v)];

    for (auto v : range(4))
        m2.vert(0, v) = old2new[m2.vert(0, v)];

    REQUIRE( m1.points.size() == 5 );

    write_by_extension("m1.geogram", m1, {{"a1", a1}, {"a2", a2}});
    write_by_extension("m2.geogram", m2, {{"a1", a1}, {"a2", a2}});
}
