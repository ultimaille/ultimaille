#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("2x2", "[svd]") {
    mat2x2 M = {{{9.5, 1.75}, {3.25, 10.45}}};
    auto [U, D, V] = svd2x2(M);
    CHECK( std::abs(D[0][0]-12.54788068808106)<1e-14 );
    CHECK( std::abs(D[1][1]-7.458430816041838)<1e-14 );
    CHECK( std::abs(D[0][1])<1e-14 );
    CHECK( std::abs(D[1][0])<1e-14 );
    CHECK( (M - U*D*V.transpose()).norm()<1e-14 );
    CHECK( std::abs(U[0]*U[1]) < 1e-14 );
    CHECK( std::abs(V[0]*V[1]) < 1e-14 );
    CHECK( std::abs(U[0].norm()-1.) < 1e-14 );
    CHECK( std::abs(U[1].norm()-1.) < 1e-14 );
    CHECK( std::abs(V[0].norm()-1.) < 1e-14 );
    CHECK( std::abs(V[1].norm()-1.) < 1e-14 );
}


