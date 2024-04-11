#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;
using namespace Linear;

TEST_CASE("test_linexpr", "[linexpr]") {
    for (int i : {0,1}) {
        LinExpr le = {{-1, 7.}, {6, -.1}, {2, -1.}, {2, 1.+1e-12}, {3, 1.}}; // test constructor
        if (i)
            le = 7 - .1*X(6) - X(2) + X(2)*(1.+1e-12) + X(3);                // test X(i) and arithmetics
        REQUIRE( le.data.size()==3 );
        CHECK((le.data[0].i == -1 && le.data[0].a == 7.));
        CHECK((le.data[1].i ==  3 && le.data[1].a == 1.));
        CHECK((le.data[2].i ==  6 && le.data[2].a == -.1));
    }
}

