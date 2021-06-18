#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

static mat<3,3> Rx(const double alpha) {
    return {{{1.,0.,0.},{0.,cos(alpha),-sin(alpha)}, {0.,sin(alpha),cos(alpha)}}};
}

static mat<3,3> Ry(const double alpha) {
    return {{{cos(alpha),0., sin(alpha)}, {0.,1.,0.},{-sin(alpha),0.,cos(alpha)}}};
}

static mat<3,3> Rz(const double alpha) {
    return {{{cos(alpha),-sin(alpha),0.}, {sin(alpha),cos(alpha),0.},{0.,0.,1.}}};
}

static double rand11() {
    return (rand()/(double)RAND_MAX)*2.-1.;
}


TEST_CASE("2x2", "[eigen decomposition]") {
    mat2x2 A = {{{4.,1.}, {1.,2.}}};
    mat2x2 evec;
    vec2 eval;
    eigendecompose_symmetric(A, eval, evec);

    vec2 eval_test = {3.+std::sqrt(2.), 3.-std::sqrt(2.)};
    mat2x2 evec_test = mat2x2{{{1./std::sqrt(4.-2.*std::sqrt(2.)), -1./std::sqrt(4.+2.*std::sqrt(2.))},{(std::sqrt(2.)-1.)/std::sqrt(4.-2.*std::sqrt(2.)), (std::sqrt(2.)+1.)/std::sqrt(4.+2.*std::sqrt(2.))}}};
    CHECK( (eval-eval_test).norm()<1e-13 );
    CHECK( (evec-evec_test).norm()<1e-13 );
}


TEST_CASE("3x3", "[eigen decomposition]") {
    mat3x3 A = {{{1.,std::sqrt(2.),2.}, {std::sqrt(2.),3.,std::sqrt(2.)},{2.,std::sqrt(2.),1.}}};
    mat3x3 evec;
    vec3 eval;
    eigendecompose_symmetric(A, eval, evec);
    CHECK( std::abs(eval[0]-5.)<1e-13 );
    CHECK((
            (std::abs(eval[1]-1.)<1e-13 && std::abs(eval[2]+1.)<1e-13) ||
            (std::abs(eval[2]-1.)<1e-13 && std::abs(eval[1]+1.)<1e-13)
         ));
}

TEST_CASE("3x3b", "[eigen decomposition]") {
    for (int test=0; test<65536; test++) {
        vec3 eval_test = {100.*rand11(), 100.*rand11(), 100.*rand11()};
        for (int i=0; i<3; i++) // sort the eigenvalues
            for (int j=i+1; j<3; j++)
                if (std::abs(eval_test[i])<std::abs(eval_test[j])) std::swap(eval_test[i], eval_test[j]);

        if (std::abs(eval_test[0])==std::abs(eval_test[1]) || std::abs(eval_test[1])==std::abs(eval_test[2])) continue; // very unlikely to happen

        mat3x3 R = Rz(2.*M_PI*rand11()) * Ry(2.*M_PI*rand11()) * Rx(2.*M_PI*rand11());
        mat3x3 A = R*mat3x3{{{eval_test[0], 0., 0.}, {0., eval_test[1], 0.}, {0., 0., eval_test[2]} }}*R.transpose();

        mat3x3 evec;
        vec3 eval;
        eigendecompose_symmetric(A, eval, evec);

        for (int i=0; i<3; i++) {
            CHECK( std::abs(eval[i]-eval_test[i])<1e-12 );
            CHECK( ((R.col(i) - evec.col(i)).norm()<1e-6 || (R.col(i) + evec.col(i)).norm()<1e-6) );
        }
    }
}

