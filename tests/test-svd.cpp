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

    CHECK( (U.invert_transpose()-U).norm() < 1e-14 );
    CHECK( (V.invert_transpose()-V).norm() < 1e-14 );
}


TEST_CASE("3x3", "[svd]") {
    mat3x3 M = {{{3,1,0}, {1,2,2}, {0,1,1}}};
    auto [U, D, V] = svd3x3(M);
//    std::cerr << U << std::endl << std::endl << D  << std::endl << std::endl << V << std::endl;
    CHECK( std::abs(D[0][0]-3.9282)<1e-4 );
    CHECK( std::abs(D[1][1]-2.3573)<1e-4 );
    CHECK( std::abs(D[2][2]-0.1079)<1e-4 );

    CHECK( std::abs(D[0][1])<1e-14 );
    CHECK( std::abs(D[0][2])<1e-14 );
    CHECK( std::abs(D[1][0])<1e-14 );
    CHECK( std::abs(D[1][2])<1e-14 );
    CHECK( std::abs(D[2][0])<1e-14 );
    CHECK( std::abs(D[2][1])<1e-14 );

    CHECK( (M - U*D*V.transpose()).norm()<1e-12 );

    CHECK( (U.invert_transpose()-U).norm() < 1e-8 );
    CHECK( (V.invert_transpose()-V).norm() < 1e-8 );
}


