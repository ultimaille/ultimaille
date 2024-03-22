#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;
using namespace Linear;

TEST_CASE("test_linexpr", "[linexpr]") {
    LinExpr le = 7 - .1*X(6) - X(2) + X(2)*(1.+1e-12) + X(3);
    REQUIRE( le.data.size()==3 );
    CHECK((le.data[0].i == -1 && le.data[0].a == 7.));
    CHECK((le.data[1].i ==  3 && le.data[1].a == 1.));
    CHECK((le.data[2].i ==  6 && le.data[2].a == -.1));
}

