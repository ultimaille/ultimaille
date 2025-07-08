#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

#include "eigen.h"
#include "svd.h"

namespace UM {

    static constexpr double precision     = 2.*std::numeric_limits<double>::epsilon();
    static constexpr double consider_null = 2.*std::numeric_limits<double>::denorm_min();

    struct Rot2D {
        double c, s;
        Rot2D operator*(const Rot2D &r) const { return {c*r.c-s*r.s, c*r.s+s*r.c}; }
        mat2x2 mat() { return {{{c,-s},{s,c}}}; }
        Rot2D transpose() const { return {c, -s}; }
    };

    void diagonalize_symmetric(Rot2D &dirot, mat2x2 &M) { // find the Jacobi rotation
        assert(std::abs(M[0][1] - M[1][0]) < 1e-14);
        if (std::abs(M[0][1])<consider_null || std::abs(M[0][0]-M[1][1])<consider_null) {
            dirot = {1.,0.};
        } else {
            double beta = (M[1][1]-M[0][0])/(2.*M[0][1]);
            double t = std::copysign(1., beta)/(std::abs(beta)+std::sqrt(beta*beta + 1.));
            dirot.c = 1./std::sqrt(t*t+1.);
            dirot.s = -dirot.c*t;
            M[0][0] -= t*M[0][1]; // M = dirot.mat().transpose()*M*dirot.mat()
            M[1][1] += t*M[0][1];
            M[0][1] = M[1][0] = 0.;
        }
    }

    void symmetrize(Rot2D &symrot, mat2x2 &M) {
        double t = M[0][0] + M[1][1];
        double d = M[1][0] - M[0][1];

        if (std::abs(t)<consider_null) {
            symrot.c = 0.;
            symrot.s = std::copysign(1., d);
        } else {
            double u = d / t;
            symrot.c = 1./std::sqrt(1.+u*u);
            symrot.s = -symrot.c*u;
        }
        M = symrot.mat() * M;
    }

    std::tuple<mat2x2,mat2x2,mat2x2> svd2x2(const mat2x2 &M) {
        mat2x2 U = mat2x2::identity();
        mat2x2 D = M;
        mat2x2 V = mat2x2::identity();

        double max_diag      = std::max(std::abs(D[1][1]), std::abs(D[0][0]));
        double max_anti_diag = std::max(std::abs(D[1][0]), std::abs(D[0][1]));

        if (max_anti_diag > std::max(consider_null, precision*max_diag)) {
            Rot2D symrot, dirot;
            symmetrize(symrot, D);
            diagonalize_symmetric(dirot, D);
            U = (symrot.transpose() * dirot).mat();
            V = dirot.mat();
        }

        for (int i : {0,1}) { // positive singular values
            if (D[i][i]<0.)
                U.set_col(i, -U.col(i));
            D[i][i] = std::abs(D[i][i]);
        }

        if (D[0][0] < D[1][1]) { // sort singular values in descending order
            std::swap(D[0][0], D[1][1]);
            std::swap(U[0][0], U[0][1]);
            std::swap(U[1][0], U[1][1]);
            std::swap(V[0][0], V[0][1]);
            std::swap(V[1][0], V[1][1]);
        }

        return {U, D, V}; // M = U * D * V.transpose()
    }

    std::tuple<mat3x3,mat3x3,mat3x3> svd3x3(const mat3x3 &M) {
        mat3x3 MtM = M.transpose() * M;
        auto [S2,V] = eigendecompose_symmetric(MtM);
        mat3x3 S = {};
        for (int i : {0,1,2}) S[i][i] = std::sqrt(S2[i]);
        mat3x3 U = M * V;
        for (int i : {0,1,2})
            for (int j : {0,1,2}) U[i][j] /= S[j][j];
        return {U, S, V};
    }
}
