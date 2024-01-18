#include "geometry.h"
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"
#include "algebra/mat.h"

namespace UM {

    // quadratures for every hex corner (https://coreform.com/papers/verdict_quality_library.pdf, p.79, 80)
    // used to compute jacobian scale on hex
    constexpr mat<8,3> QH[9] = { 
        {{
            // A_0 = (L_0, L_3, L_4)
            {-1,-1,-1},
            {1,0,0},
            {0,0,0},
            {0,1,0},
            {0,0,1},
            {0,0,0},
            {0,0,0},
            {0,0,0}
        }},
        {{
            // A_1 = (L_1, -L_0, L_5)
            {0,1,0},
            {-1,-1,-1},
            {1,0,0},
            {0,0,0},
            {0,0,0},
            {0,0,1},
            {0,0,0},
            {0,0,0}
        }},
        {{
            // A_2 = (L_2, -L_1, L_6)
            {0,0,0},
            {0,1,0},
            {-1,-1,-1},
            {1,0,0},
            {0,0,0},
            {0,0,0},
            {0,0,1},
            {0,0,0}
        }},
        {{
            // A_3 = (-L_3, -L_2, L_7)
            {1,0,0},
            {0,0,0},
            {0,1,0},
            {-1,-1,-1},
            {0,0,0},
            {0,0,0},
            {0,0,0},
            {0,0,1}
        }},
        {{
            // A_4 = (L_11, L_8, -L_4)
            {0,0,1},
            {0,0,0},
            {0,0,0},
            {0,0,0},
            {-1,-1,-1},
            {0,1,0},
            {0,0,0},
            {1,0,0}
        }},
        {{
            // A_5 = (-L_8, L_9, -L_5)
            {0,0,0},
            {0,0,1},
            {0,0,0},
            {0,0,0},
            {1,0,0},
            {-1,-1,-1},
            {0,1,0},
            {0,0,0}
        }},
        {{
            // A_6 = (-L_9, L_10, -L_6)
            {0,0,0},
            {0,0,0},
            {0,0,1},
            {0,0,0},
            {0,0,0},
            {1,0,0},
            {-1,-1,-1},
            {0,1,0}
        }},
        {{
            // A_7 = (-L_10, -L_11, -L_7)
            {0,0,0},
            {0,0,0},
            {0,0,0},
            {0,0,1},
            {0,1,0},
            {0,0,0},
            {1,0,0},
            {-1,-1,-1}
        }},
        {{
            // A_8 = (X_1, X_2, X_3) = {(P_1-P_0)+(P_2-P_3)+(P_5-P_4)+(P_6-P_7)}
            {-1,-1,-1},
            {1,-1,-1},
            {1,1,-1},
            {-1,1,-1},
            {-1,-1,1},
            {1,-1,1},
            {1,1,1},
            {-1,1,1}
        }}
    };    

    // quadratures for every hex corner (counter-clock wise)
    // used to compute jacobian scale on quad
    constexpr mat<4,2> QQ[4] = { 
        {{
            // A_0 = (L_0, L_3)
            {-1,-1},
            {1,0},
            {0,0},
            {0,1},
        }},
        {{
            // A_1 = (L_0, L_1)
            {-1,0},
            {1,-1},
            {0,1},
            {0,0},
        }},
        {{
            // A_2 = (-L_2, L_1)
            {0,0},
            {0,-1},
            {1,1},
            {-1,0},
        }},
        {{
            // A_3 = (-L_2, -L_3)
            {0,-1},
            {0,0},
            {1,0},
            {-1,1},
        }}
    };

    Triangle2 Triangle3::project() const {
        const vec3 &A = v[0];
        const vec3 &B = v[1];
        const vec3 &C = v[2];

        vec3 X = (B - A).normalized(); // construct an orthonormal 3d basis
        vec3 Z = cross(X, C - A).normalized();
        vec3 Y = cross(Z, X);

        
        const vec2 z0 = vec2(0,0); // project the triangle to the 2d basis (X,Y)
        const vec2 z1 = vec2((B - A).norm(), 0);
        const vec2 z2 = vec2((C - A)*X, (C - A)*Y);

        return Triangle2{{z0, z1, z2}};
    }

    double Quad3::unsigned_area() const {
        return UM::unsigned_area(v, 4);
    }

    vec3 Quad3::normal() const {
        return UM::normal(v, 4);
    }

    vec3 Quad3::bary_verts() const {
        return UM::bary_verts(v, 4);
    }

    double Quad3::jacobian(int c) const {

        mat<2,4> A = {{
            {v[0].x, v[1].x, v[2].x, v[3].x},
            {v[0].y, v[1].y, v[2].y, v[3].y}
        }};

        const mat<4,2> &B = QQ[c];

        mat2x2 J = A * B;
        // scale
        const vec2& l0 = J.col(0);
        const vec2& l1 = J.col(1);

        // Check L_min² <= DBL_MIN => q = DBL_MAX, with DBL_MIN = 1e-30 and DBL_MAX = 1e+30
        const double l_min = std::min(l0.norm2(), l1.norm2());

        if (l_min <= 1e-30)
            return 1e30;

        double d = std::sqrt(l0.norm2() * l1.norm2());

        double scaled_jacobian = J.det() / d;

        if (scaled_jacobian > 0)
            return std::min(scaled_jacobian, 1e30);
        else 
            return std::max(scaled_jacobian, -(1e30));
    }

    double Quad3::scaled_jacobian() const {
        double min_j = std::numeric_limits<double>::max();
        for (int c = 0; c < 4; c++) {
            double j = jacobian(c);
            if (j < min_j)
                min_j = j;
        }

        return min_j;
    }

    vec3 Poly3::bary_verts() const {
        return UM::bary_verts(v.data(), static_cast<const int>(v.size()));
    }

    vec3 Poly3::normal() const {
        return UM::normal(v.data(), static_cast<const int>(v.size()));
    }

    double Poly3::unsigned_area() const {
        return UM::unsigned_area(v.data(), static_cast<const int>(v.size()));
    }

    vec3 Tetrahedron3::bary_verts() const {
        return UM::bary_verts(v, 4);
    }

    vec3 Hexahedron3::bary_verts() const {
        return UM::bary_verts(v, 8);
    }

    double Hexahedron3::volume() const {
        const vec3 bary = bary_verts();
        double vol = 0;

        // Indexes of hex faces
        // Vertices numbering convention (from geogram) is very IMPORTANT (and respected)
        // as well as counter-clock wise convention for facet orientation
        // TODO use volume_reference.h instead
        int indexes[] = {
            0,1,5,4, // front
            2,6,7,3, // back
            2,0,4,6, // left
            1,3,7,5, // right
            4,5,7,6, // top
            0,2,3,1 // bottom
        };

        for (int f=0; f<6; f++) {
            for (int i=0; i<4; i++) {
                const int offset = 4*f;
                vol += tet_volume(
                        bary,
                        v[indexes[i + offset]],
                        v[indexes[(i+1)%4 + offset]],
                        v[indexes[(i+2)%4 + offset]]
                        )*.5;
            }
        }
        return vol;
    }

    // // evaluate the Jacobian matrix at the given corner, return the Jacobian determinant
    // double Hexahedron3::jacobian(int c) const {
    //     double J[3][3];

    //     double A[3][8] = {
    //         {v[0].x, v[1].x, v[2].x, v[3].x, v[4].x, v[5].x, v[6].x, v[7].x},
    //         {v[0].y, v[1].y, v[2].y, v[3].y, v[4].y, v[5].y, v[6].y, v[7].y},
    //         {v[0].z, v[1].z, v[2].z, v[3].z, v[4].z, v[5].z, v[6].z, v[7].z}
    //     };

    //     const double (&B)[8][3] = Q[c];

    //     // TODO can compute using mat
    //     //J = A * B;
    //     for (int j=0; j<3; j++) // J = A*B
    //         for (int i=0; i<3; i++) {
    //             J[j][i] = 0;
    //             for (int k=0; k<8; k++)
    //                 J[j][i] += A[j][k]*B[k][i];
    //         }



    //     // 3x3 det (TODO can compute using v1.(v2Xv3))
    //     return (J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1]) -
    //         (J[0][2]*J[1][1]*J[2][0] + J[0][1]*J[1][0]*J[2][2] + J[0][0]*J[1][2]*J[2][1]);
    // }

    double Hexahedron3::jacobian(int c) const {

        // Convert verdict volume reference to ultimaille volume reference (https://coreform.com/papers/verdict_quality_library.pdf, p.79)
        // Note: Only need to invert points at index 2,3 and 6,7
        mat<3,8> A = {{
            {v[0].x, v[1].x, v[3].x, v[2].x, v[4].x, v[5].x, v[7].x, v[6].x},
            {v[0].y, v[1].y, v[3].y, v[2].y, v[4].y, v[5].y, v[7].y, v[6].y},
            {v[0].z, v[1].z, v[3].z, v[2].z, v[4].z, v[5].z, v[7].z, v[6].z}
        }};

        const mat<8,3> &B = QH[c];

        mat3x3 J = A * B;
        // scale
        const vec3& l0 = J.col(0);
        const vec3& l1 = J.col(1);
        const vec3& l2 = J.col(2);

        // Check L_min² <= DBL_MIN => q = DBL_MAX, with DBL_MIN = 1e-30 and DBL_MAX = 1e+30
        const double l_min = std::min(std::min(l0.norm2(), l1.norm2()), l2.norm2());

        if (l_min <= 1e-30)
            return 1e30;

        double d = std::sqrt(l0.norm2() * l1.norm2() * l2.norm2());

        double scaled_jacobian = J.det() / d;

        if (scaled_jacobian > 0)
            return std::min(scaled_jacobian, 1e30);
        else 
            return std::max(scaled_jacobian, -(1e30));
    }

    double Hexahedron3::scaled_jacobian() const {
        double min_j = std::numeric_limits<double>::max();
        for (int c = 0; c < 9; c++) {
            double j = jacobian(c);
            if (j < min_j)
                min_j = j;
        }

        return min_j;
    }

    vec3 Pyramid3::bary_verts() const {
        return UM::bary_verts(v, 5);
    }

    double Pyramid3::volume() const {

        // Indexes of pyramid base
        // Vertices numbering convention (from geogram) is very IMPORTANT
        // as well as counter-clock wise convention for facet orientation
        int indexes[] = {3,2,1,0};

        const vec3 apx = apex();
        double vol = 0;

        // Compute base to apex volume
        for (int i=0; i<4; i++) {
            vol += tet_volume(
                    apx,
                    v[indexes[i]],
                    v[indexes[(i+1)%4]],
                    v[indexes[(i+2)%4]]
                    )*.5;
        }

        return vol;
    }

}