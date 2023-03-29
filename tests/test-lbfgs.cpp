#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

static const double xtol = 1e-3;

TEST_CASE("2d quadratic", "[L-BFGS]") {
    std::vector<double> sol = {10, 10};

    const auto func = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = (x[0] - 7)*(x[0] - 7) +
            (x[1] - 1)*(x[1] - 1);
        g[0] = 2*(x[0] - 7);
        g[1] = 2*(x[1] - 1);
    };

    STLBFGS::Optimizer opt(func);
    opt.run(sol);

    CHECK(std::abs(sol[0]-7)<xtol);
    CHECK(std::abs(sol[1]-1)<xtol);
}

