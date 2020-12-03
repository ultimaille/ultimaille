#include <iostream>

#include <ultimaille/all.h>

using namespace UM;

int main() {
    std::vector<double> sol(2);
    sol[0] = 10;
    sol[1] = 10;
    auto func = [](std::vector<double>& x, double& f, std::vector<double>& g) {
        f = (x[0] - 7) * (x[0] - 7) + (x[1] - 1) * (x[1] - 1);
        g[0] = 2 * x[0] - 14;
        g[1] = 2 * x[1] - 2;
    };

#if 0
    hlbfgs_optimizer::optimize(func, sol);
#else
    hlbfgs_optimizer opt(func);
    opt.set_verbose(true);
    opt.optimize(sol);
#endif

    std::cerr << sol[0] << " " << sol[1] << std::endl;
    return 0;
}

