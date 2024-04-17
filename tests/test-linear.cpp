#include <catch2/catch_test_macros.hpp>
#include <ultimaille/all.h>

using namespace UM;
using namespace Linear;

TEST_CASE("test_linexpr", "[linexpr]") {
    for (int i : {0,1}) {
        LinExpr le = {{-1, 7.}, {6, -.1}, {2, -1.}, {2, 1.+1e-12}, {3, 1.}}; // test constructor
        if (i)
            le = 7 - .1*X(6) - X(2) + X(2)*(1.+1e-12) + X(3);                // test X(i) and arithmetics
        REQUIRE( le.data.size()==3 );
        CHECK((le.data[0].index == -1 && le.data[0].value == 7.));
        CHECK((le.data[1].index ==  3 && le.data[1].value == 1.));
        CHECK((le.data[2].index ==  6 && le.data[2].value == -.1));
    }
}

