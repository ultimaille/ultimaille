#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

static const std::string points1_str =
R"(-4.137 -2.376 7.58324
-3.78244 -2.93098 9.43324
-3.67289 -2.63485 9.35324
-3.68462 -2.52142 8.85324
-3.47944 -3.18353 6.65324
-1.4844 -3.62345 7.12324
-1.35024 -3.30096 7.12324
-1.13755 -3.6096 6.49324
-0.999408 -3.63234 7.10324
-1.12841 -3.24627 7.10324
-1.30205 -3.56226 7.07324
-1.41606 -3.26984 6.72324
-1.54896 -3.39997 7.10324
)";

TEST_CASE("PointSet XYZ IO test", "[XYZ]") {
    static const std::string filename[2] = { "ultimaille-test-pointset-in.xyz", "ultimaille-test-pointset-out.xyz" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << points1_str;
    ofs.close();

    PointSet m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].size()==13 );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
}


