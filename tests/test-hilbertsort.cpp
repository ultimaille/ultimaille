#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

static double rand01() {
    return (rand()/(double)RAND_MAX);
}

TEST_CASE("Test max dist", "[Hilbert sort]") {
    int n = 1000000;
    PointSet pts;

    for ([[maybe_unused]] int i : range(n))
        pts.push_back({rand01(), rand01(), rand01()});

    Permutation perm(n);

    HilbertSort hs(*pts.data);
    hs.apply(perm.ind);

    perm.apply(*pts.data);

    double maxdist = -std::numeric_limits<double>::max();
    for (int i=0; i<n-1; i++)
        maxdist = std::max(maxdist, (pts[i] - pts[i+1]).norm());

    CHECK( 1./std::cbrt(n) * 7 > maxdist ); // n^(1/3) samples per cube edge, ideally n^(-1/3) max dist, take 7 times that to be sure
}

