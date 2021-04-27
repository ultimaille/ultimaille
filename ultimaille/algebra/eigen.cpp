#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "eigen.h"

namespace UM {

    // Let A be a 2x2 real symmetric matrix, then A has all real eigenvalues and a full set of orthogonal eigenvectors.
    // We want to compute A = R D R^T, where R = [v1 v2] is an orthogonal matrix whose columns are the eigenvectors,
    // and D is the diagonal matrix built with the eigenvalues.
    // A = [a11 a12]   R = [cos t -sin t]
    //     [a12 a22]       [sin t  cos t]
    // Write down D = R^T A R, and impose the diagonality constraint on it (equate to zero the off-diagonal element).
    // It allows us to compute (cos 2t, sin 2t). Then deduce (cos t, sin t).

    void eigendecompose_symmetric(const double a11, const double a12, const double a22, vec2 &eval, mat2x2 &evec) {
        double cos2t = .5*(a11 - a22), sin2t = a12;  // (cos 2t, sin 2t) vector, needs normalization to be actual sine and cosine
        // (cos 2t, sin 2t) is defined up to a sign: choose the sign in a way that cos(2t)<=0, it helps numerical stability for computing sin(t)
        double maxabs = std::max(std::abs(cos2t), std::abs(sin2t));
        if (maxabs > 0.) { // take care of overflows/underflows
            cos2t /= maxabs;  // in [-1,1]
            sin2t /= maxabs;  // in [-1,1]
            double length = std::sqrt(cos2t*cos2t + sin2t*sin2t);
            cos2t /= length;
            sin2t /= length;
            if (cos2t > 0.) {
                cos2t = -cos2t;
                sin2t = -sin2t;
            }
        } else {
            cos2t = -1.;
            sin2t =  0.;
        }

        // now compute (cos t, sin t)
        double sint = std::sqrt(.5*(1. - cos2t));  // cos(2t)<=0, therefore sin(t) >= 1/sqrt(2)
        double cost = .5*sin2t/sint;                 // no division by zero here

        eval[0] = cost*cost*a11 + sin2t*a12 + sint*sint*a22;
        eval[1] = cost*cost*a22 - sin2t*a12 + sint*sint*a11;

        if (std::abs(eval[0])>std::abs(eval[1])) {
            evec = { {{cost, -sint}, {sint, cost}} }; // vectors are the columns
        } else {
            std::swap(eval[0], eval[1]);
            evec = { {{sint, cost}, {-cost, sint}} }; // permute the vectors and negate one of them to keep the basis right-handed
        }
    }

    void eigendecompose_symmetric(const mat2x2 &A, vec2 &eval, mat2x2 &evec) {
        assert(std::abs(A[0][1]-A[1][0])<1e-13);
        eigendecompose_symmetric(A[0][0], A[0][1], A[1][1], eval, evec);
    }


    // Given a symmetric matrix, compute a quaternion q such that its corresponding matrix Q can be used to diagonalize A
    // i.e. A = Q D Q^T, where D is the diagonal matrix built with the eigenvalues: D = Q^T A Q
    Quaternion eigenvectors_symmetric(const mat3x3 &A) {
        for (int i=0; i<2; i++) // assert the symmetry
            for (int j=i+1; j<3; j++)
                assert(std::abs(A[i][j] - A[j][i])<1e-13);

        Quaternion Q = {{0.,0.,0.}, 1.}; // zero rotation
        for (int step=0; step<64; step++) { // certainly won't need that many steps
            mat3x3 M = Q.rotation_matrix();
            mat3x3 D = M.transpose() * A * M;

            int p = 1, q = 2; // find the largest magnitude off-diagonal element of D
            for (int j=1; j<3; j++)
                if (std::abs(D[p][q]) < std::abs(D[0][j])) {
                    p = 0;
                    q = j;
                }

            if (D[p][q]==0.) break; // the matrix is already diagonal

            double cot2phi = (D[q][q]-D[p][p])/(2.*D[p][q]);
            double sgn = cot2phi > 0. ? 1. : -1.;
            cot2phi = std::abs(cot2phi);

            double tanphi = sgn/(cot2phi + (cot2phi>1e20 ? cot2phi : std::sqrt(cot2phi*cot2phi + 1.)));
            double cosphi = 1./std::sqrt(tanphi*tanphi + 1.);

            if (cosphi==1.) break;  // no room for improvement - reached machine precision.

            Quaternion J = {{0.,0.,0.}, 1.}; // Jacobi rotation for this iteration
            int axis = 3-p-q;        // quaternion axis for the Jacobi rotation
            if (1==axis) sgn = -sgn; // since quaternion rotations naturally follow the right-hand rule, change the sign to be sure that Jpq = -sin(phi)
            J[axis] = -sgn*std::sqrt((1.-cosphi)/2.); // use half-angle formula to obtain sin(phi/2) from cos(phi)
            J.w = std::sqrt((1.+cosphi)/2.);          // and cos(phi/2) from cos(phi)

            if (J.w==1.) break; // reached limits of floating point precision
            Q = Q*J;
        }
        return Q.normalize();
    }

    void eigendecompose_symmetric(const mat3x3 &A, vec3 &eval, mat3x3 &evec) {
        Quaternion Q = eigenvectors_symmetric(A);
        mat3x3 M = (Q.rotation_matrix()).transpose(); // to ease access to the eigenvectors
        mat3x3 D = M * A * M.transpose();
        eval = { D[0][0], D[1][1], D[2][2] };

        if (std::abs(eval[0])<std::abs(eval[1])) { // sort the eigevnalues by the magnitude
            std::swap(eval[0], eval[1]);
            std::swap(M[0], M[1]);
        }
        if (std::abs(eval[0])<std::abs(eval[2])) {
            std::swap(eval[0], eval[2]);
            std::swap(M[0], M[2]);
        }
        if (std::abs(eval[1])<std::abs(eval[2])) {
            std::swap(eval[1], eval[2]);
            std::swap(M[1], M[2]);
        }
        if (cross(M[0], M[1])*M[2]<0) { // be sure that the eigenvectors form a right-hand basis
            M[2] = -M[2];
        }
        evec = M.transpose();
    }

}

