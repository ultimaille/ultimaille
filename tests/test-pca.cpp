#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch_test_macros.hpp>
//#include <iostream>
//#include <iomanip>
//#include <limits>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("3x3 PCA", "PointSet") {
    PointSet p;
    *p.data = { {90., 60., 90.}, {90., 90., 30.}, {60., 60., 60.}, {60., 60., 90.}, {30., 30., 30.} };
    PointSetCovariance cov(*p.data);

    mat3x3 C = {{{504,360,180}, {360,360,0}, {180,0,720}}};

    CHECK ( (cov.cov - C).norm() < 1e-10 );
    CHECK ( (cov.center - vec3{66,60,60}).norm() < 1e-10 );

    auto [axes, eval, center] = Inspect(p).principal_axes();

    CHECK ( (center - vec3{66,60,60}).norm() < 1e-10 );
    CHECK( (eval - vec3{910.069953, 629.110387, 44.819660}).norm() < 1e-5 );
    CHECK( (axes - mat3x3{{{0.655802, 0.429198, 0.621058}, {-0.385999, -0.516366, 0.764441}, {0.648790, -0.741050, -0.172964}}}.transpose()).norm() < 1e-5 );
}

TEST_CASE("3x3 covariance", "PointSet") {
    std::vector<vec3> A(10);
    std::vector<vec3> B(20);
    for (auto& P : A) for (int d : {0,1,2}) P[d] = (d+1) * (rand() % 100);
    for (auto& P : B) for (int d : {0,1,2}) P[d] = (d+1) * (rand() % 100);
    std::vector<vec3> AB;
    for (auto& P : A) AB.push_back(P);
    for (auto& P : B) AB.push_back(P);

    PointSetCovariance mA(A);
    PointSetCovariance mB(B);
    PointSetCovariance mAB(AB);
    PointSetCovariance mAunionB = mA + mB;
//  std::cerr << std::setprecision(std::numeric_limits<double>::max_digits10);
//  std::cerr << mAB.cov << std::endl << mAunionB.cov << std::endl << (mAB.cov - mAunionB.cov).norm() << std::endl;
    CHECK( (mAB.cov - mAunionB.cov).norm() < 1e-9 );
}

