#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

// vt per vertex, vn present
static const std::string obj1_str =
R"(
v -1 0 -1
v  1 0 -1
v  1 0  1
v -1 0  1

vt 0 0
vt 1 0
vt 1 1
vt 0 1

vn 0 1 0

f 3/3/1 2/2/1 1/1/1
f 4/4/1 3/3/1 1/1/1
)";

// vt per vertex, vn absent
static const std::string obj2_str =
R"(
v -1 0 -1
v  1 0 -1
v  1 0  1
v -1 0  1

vt 0 0
vt 1 0
vt 1 1
vt 0 1

f 3/3 2/2 1/1
f 4/4 3/3 1/1
)";

// vt per corner
static const std::string obj3_str =
R"(
v -1 0 -1
v  1 0 -1
v  1 0  1
v -1 0  1

vt 0 0
vt 1 0
vt 1 1
vt 0 1

f 3/1 2/2 1/3
f 4/1 3/3 1/4
)";


// vt per corner, vn present, 1 quad + 2 triangles
static const std::string obj4_str =
R"(
v -1 0 -1
v  1 0 -1
v  1 0  1
v -1 0  1

v -1 1 -1
v  1 1 -1
v  1 1  1
v -1 1  1

vt 0 0
vt 1 0
vt 1 1
vt 0 1

vn 0 1 0

f 3/3/1 2/2/1 1/1/1
f 4/4/1 3/3/1 1/1/1
f 5/1/1 6/2/1 7/3/1 8/4/1
)";


TEST_CASE("tex_coord per vertex IO test#1", "[OBJ]") {
    static const std::string filename[2] = { "ultimaille-test-tex_coord1-in.obj", "ultimaille-test-tex_coord1-out.obj" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << obj1_str;
    ofs.close();

    Triangles m[2] = {};
    for (int i : range(2)) {
        SurfaceAttributes attr = read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].nfacets()==2 );

        int cnt = 0;
        for (auto &pair : attr.points)
            cnt += (pair.first=="tex_coord");
        CHECK( cnt==1 );

        for (int f : range(2))
            CHECK( std::abs(m[i].util.unsigned_area(f)-2.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0], attr);
    }
}

TEST_CASE("tex_coord per vertex IO test#2", "[OBJ]") {
    static const std::string filename[2] = { "ultimaille-test-tex_coord2-in.obj", "ultimaille-test-tex_coord2-out.obj" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << obj2_str;
    ofs.close();

    Triangles m[2] = {};
    for (int i : range(2)) {
        SurfaceAttributes attr = read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].nfacets()==2 );

        int cnt = 0;
        for (auto &pair : attr.points)
            cnt += (pair.first=="tex_coord");
        CHECK( cnt==1 );

        for (int f : range(2))
            CHECK( std::abs(m[i].util.unsigned_area(f)-2.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0], attr);
    }
}

TEST_CASE("tex_coord per corner IO test", "[OBJ]") {
    static const std::string filename[2] = { "ultimaille-test-tex_coord3-in.obj", "ultimaille-test-tex_coord3-out.obj" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << obj3_str;
    ofs.close();

    Triangles m[2] = {};
    for (int i : range(2)) {
        SurfaceAttributes attr = read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].nfacets()==2 );

        int cnt = 0;
        for (auto &pair : attr.corners)
            cnt += (pair.first=="tex_coord");
        CHECK( cnt==1 );

        for (int f : range(2))
            CHECK( std::abs(m[i].util.unsigned_area(f)-2.)<ftol );

        if (!i)
            write_by_extension(filename[1], m[0], attr);
    }
}

TEST_CASE("Polygons IO test", "[OBJ]") {
    static const std::string filename[2] = { "ultimaille-test-tex_coord4-in.obj", "ultimaille-test-tex_coord4-out.obj" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << obj4_str;
    ofs.close();

    Polygons m[2] = {};
    for (int i : range(2)) {
        SurfaceAttributes attr = read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==8 );
        REQUIRE( m[i].nfacets()==3 );

        int cnt = 0;
        for (auto &pair : attr.corners)
            cnt += (pair.first=="tex_coord");
        CHECK( cnt==1 );

//      for (int f : range(2))
//          CHECK( std::abs(m[i].util.unsigned_area(f)-2.)<ftol );
//      CHECK( std::abs(m[i].util.unsigned_area(2)-4.)<ftol );

        if (!i)
            write_by_extension(filename[1], m[0], attr);
    }
}

