#include <catch2/catch_test_macros.hpp>
#include <ultimaille/all.h>

using namespace UM;
using namespace Linear;

TEST_CASE("test_linexpr", "[LinExpr]") {
    for (int i : {0,1}) {
        LinExpr le = {{-1, 7.}, {6, -.1}, {2, -1.}, {2, 1.+1e-12}, {3, 1.}}; // test constructor
        if (i)
            le = 7 - .1*X(6) - X(2) + X(2)*(1.+1e-12) + X(3);                // test X(i) and arithmetics
        REQUIRE( le.size()==3 );
        CHECK((le[0].index == -1 && le[0].value == 7.));
        CHECK((le[1].index ==  3 && le[1].value == 1.));
        CHECK((le[2].index ==  6 && le[2].value == -.1));
    }
}

TEST_CASE("SparseVector operations", "[SparseVector]") {
    SECTION("Compact") {
        SparseVector vec = {{1, 3.1}, {3, 2.}, {5, -1e-14}, {6, 10.}, {3, 1.5}};
        REQUIRE(vec.size() == 3);
        vec *= 3.14;
        REQUIRE((vec[0].index==1 && std::abs(vec[0].value - 3.14*3.1)<1e-10));
        REQUIRE((vec[1].index==3 && std::abs(vec[1].value - 3.14*3.5)<1e-10));
        REQUIRE((vec[2].index==6 && std::abs(vec[2].value - 3.14*10.)<1e-10));
    }
}

TEST_CASE("LOLMatrix methods", "[LOLMatrix]") {
    LOLMatrix m = {{
        {{0, 1.0}, {2, 2.0}, {3, 3.0}},
        {{1, 1.0}, {7, 2.0}},
        {{0, 1.0}, {2, 2.0}, {2, 3.0}, {3, 4.0}}
    }};

    REQUIRE(m.nrows() == 3);
    REQUIRE(m.count_columns() == 8);
    REQUIRE(m.count_nnz() == 8);

    m.drop_zero_columns();

    CRSMatrix crs = m.to_crs();
    std::vector<double> crsv = {1., 2., 3., 1., 2., 1., 5., 4.};
    std::vector<int>    crsi = {0, 2, 3, 1, 4, 0, 2, 3};

    REQUIRE( (crs.mat.size() == 8 && crs.nnz() == 8) );
    REQUIRE( crs.nrows() == 3 );
    for (int i=0; i<8; i++) {
        REQUIRE( crs.mat[i].value == crsv[i] );
        REQUIRE( crs.mat[i].index == crsi[i] );
    }
}


