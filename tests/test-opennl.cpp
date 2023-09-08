#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

static const double xtol = 1e-3;

TEST_CASE("test_simple_linear_solve", "[opennl]") {
// solve following linear system:
// 1.0*x0 + 2.0*x1 = 5.0
// 3.0*x0 + 4.0*x1 = 6.0

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, 2);
    nlSolverParameteri(NL_SOLVER, NL_SOLVER_DEFAULT);

    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);
    nlBegin(NL_ROW);
    nlCoefficient(0, 1.);
    nlCoefficient(1, 2.);
    nlRightHandSide(5.);
    nlEnd(NL_ROW);
    nlBegin(NL_ROW);
    nlCoefficient(0, 3.);
    nlCoefficient(1, 4.);
    nlRightHandSide(6.);
    nlEnd(NL_ROW);
    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);

    nlSolve();

    CHECK(std::abs(1.*nlGetVariable(0) + 2.*nlGetVariable(1) - 5.)<xtol);
    CHECK(std::abs(3.*nlGetVariable(0) + 4.*nlGetVariable(1) - 6.)<xtol);

    nlDeleteContext(nlGetCurrent());
}

