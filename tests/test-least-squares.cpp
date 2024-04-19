#include <catch2/catch_test_macros.hpp>
#include <ultimaille/all.h>

using namespace UM;

static const double xtol = 1e-3;

TEST_CASE("unconstrained", "[LeastSquares]") {
    LeastSquares ls(3);
    ls.fix(0, 17.);
    ls.add_to_energy(Linear::X(1) - 7.);
    ls.add_to_energy(Linear::X(2) - 3.);
    ls.solve();

    CHECK(std::abs(ls.value(0)-17.)<xtol);
    CHECK(std::abs(ls.value(1)- 7.)<xtol);
    CHECK(std::abs(ls.value(2)- 3.)<xtol);
}

TEST_CASE("constrained", "[ConstrainedLeastSquares]") {
    using namespace Linear;
    ConstrainedLeastSquares cls(4);
    cls.add_to_constraints(X(0)+X(1)+X(2)-X(3)-1);
    cls.add_to_constraints(X(0)-X(1)+X(2)+X(3)-3);
    cls.add_to_constraints(X(0)+X(1)-X(2)+X(3)+1);

    cls.add_to_energy(X(0)+X(1)+X(2)+X(3) - 2);
    cls.add_to_energy(X(0)+3*X(1)+X(2)+X(3) - 1);
    cls.add_to_energy(X(0)-X(1)+3*X(2)+X(3) - 6);
    cls.add_to_energy(X(0)+X(1)+X(2)+3*X(3) - 3);
    cls.add_to_energy(X(0)+X(1)+X(2)-X(3) - 1);

    cls.solve();

    CHECK(std::abs(cls.value(0)-0.5)<xtol);
    CHECK(std::abs(cls.value(1)+0.5)<xtol);
    CHECK(std::abs(cls.value(2)-1.5)<xtol);
    CHECK(std::abs(cls.value(3)-0.5)<xtol);
}


