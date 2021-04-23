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

static const double ftol = 1e-13;

TEST_CASE("Quaternion", "[quaternions]") {
    const double alpha = M_PI/6.;
    const Quaternion q1 = {vec3{0.,0.,1.}*sin(alpha/2.), cos(alpha/2.)};
    const mat<3,3> M1 = Rz(alpha);

    { // test quaternion->rotation matrix
        mat<3,3> R = q1.rotation_matrix();
        CHECK( (M1-R).norm()<ftol );
    }

    const Quaternion q2 = {{0.592817248117098, 0.083109566226999, 0.597780725760344}, 0.533215448243828};
    const mat<3,3> M2 = {{{0.271502007821992, -0.538954266589856, 0.797380058863537}, {0.736029403961437, -0.417548171511386, -0.532836035729257}, {0.620118840427219, 0.731561222996468, 0.283321020672862}}};

    { // test quaternion->rotation matrix
        mat<3,3> R = q2.rotation_matrix();
        CHECK((M2-R).norm()<ftol);
    }

    { // test multiplication of quaternions
        Quaternion q3 = q2*q1;
        mat<3,3> R = q3.rotation_matrix();
        CHECK((M2*M1-R).norm()<ftol);
    }

    { // test quaternion->Euler angles
        vec3 ea = q2.euler_angles();
        mat<3,3> M = Rz(ea.z)*Ry(ea.y)*Rx(ea.x);
        mat<3,3> R = q2.rotation_matrix();
        CHECK((M-R).norm()<ftol);
    }
}

