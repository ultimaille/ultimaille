#include "geometry.h"
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"
#include "algebra/mat.h"
#include "volume_reference.h"
#include "volume.h"

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

    // quadratures for every quad corner (counter-clock wise)
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

	vec3 Triangle2::bary_coords(vec2 G) const {
		vec3 result;
		double sum = 0;

        for (int i = 0; i < 3; i++) {
			const vec2 &A = v[(i + 1) % 3];
			const vec2 &B = v[(i + 2) % 3];

			result[i] = Triangle2{{A, B, G}}.signed_area();
			sum += result[i];
		}

		return result / sum;
	}

    vec3 Triangle3::bary_coords(vec3 G) const {
        vec3 result;
        double sum = 0;
        const vec3 n = normal();
        
        for (int i = 0; i < 3; i++) {
            const vec3 &A = v[(i + 1) % 3];
            const vec3 &B = v[(i + 2) % 3];

            result[i] = Triangle3{{A, B, G}}.cross_product() * n;
            sum += result[i];
        }

        return result / sum;
    }

    vec4 Tetrahedron3::bary_coords(vec3 G) const {
        vec4 result;
        double sum = 0;
        
        const double vol = volume();

        for (int i = 0; i < 4; i++) {
            const vec3 &A = v[i % 4];
            const vec3 &B = v[(i + 1) % 4];
            const vec3 &C = v[(i + 2) % 4];

            result[i] = std::abs(Tetrahedron3{{A, B, C, G}}.volume() / vol);
            sum += result[i];
        }

        return result / sum;
    }

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

    Triangle2 Triangle2::dilate(double scale) const {
        vec2 G = bary_verts();
        return {{G + scale * (v[0] - G), G + scale * (v[1] - G), G + scale * (v[2] - G)}};
    }

    Triangle3 Triangle3::dilate(double scale) const {
        vec3 G = bary_verts();
        return {{G + scale * (v[0] - G), G + scale * (v[1] - G), G + scale * (v[2] - G)}};
    }
    
    mat3x3 Triangle3::tangent_basis() const {
        mat3x3 res = {v[1] - v[0], v[2] - v[0], vec3(0,0,0)};
        
        for (int d = 0; d < 2; d++)
            res[d].normalize();
        
        res[2] = cross(res[0], res[1]);
        res[2].normalize();
        res[1] = cross(res[2], res[0]);

        return res;
    }

    mat3x3 Triangle3::tangent_basis(vec3 first_axis) const {
        mat3x3 res = {v[1] - v[0], v[2] - v[0]};
        
        for (int d = 0; d < 2; d++)
            res[d].normalize();
        
        res[2] = cross(res[0], res[1]);
        res[2].normalize();
        res[0] = first_axis;
        res[0].normalize();
        res[1] = cross(res[2], res[0]);
        return res;
    }

    mat3x3 Triangle3::grad_operator() const {
        mat3x3 accum;
        vec3 tr_normal = normal();
        
        for (int i = 0; i < 3; i++) {
            const vec3 &P = v[i];
            const vec3 &A = v[(i + 1) % 3];
            const vec3 &B = v[(i + 2) % 3];
            vec3 n = cross(B - A, tr_normal);
            n.normalize();
            double scale = (n * (P - A));
            accum[i] = n / scale;
        }

        return accum.transpose();
    }

    vec3 Triangle3::grad(vec3 u) const {
        return grad_operator() * u;
    }

	mat<2,3> Triangle2::grad_operator() const {
		mat<3,2>  accum;

        for (int i = 0; i < 3; i++) {
			const vec2 &P = v[i];
			const vec2 &A = v[(i + 1) % 3];
			const vec2 &B = v[(i + 2) % 3];

			vec2 n((B - A)[1], -(B - A)[0]);
			n.normalize();
			double scale = (n * (P - A));
			accum[i] = n / scale;
		}
		return accum.transpose();
	}

	vec2 Triangle2::grad(vec3 u) const {
		return grad_operator() * u;
	}

    mat<3,4> Tetrahedron3::grad_operator() const {
        mat<4,3>  accum ;

        for (int i = 0; i < 3; i++) {
            const vec3 &P = v[i];
            const vec3 &A = v[(i + 1) % 4];
            const vec3 &B = v[(i + 2) % 4];
            const vec3 &C = v[(i + 3) % 4];
            vec3 n = cross(B - A, C - A);
            n.normalize();
            double scale = (n * (P - A));
            accum[i] = n / scale;
        }

        return accum.transpose();
    }

    vec3 Tetrahedron3::grad(vec4 u) const {
        return grad_operator() * u;
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

    double Quad2::jacobian(int c) const {

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

    double Quad2::scaled_jacobian() const {
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
        // Vertices numbering convention is very IMPORTANT (and respected in volume_reference.h)
        // as well as counter-clock wise convention for facet orientation
        const auto &indexes = reference_cells[Volume::CELL_TYPE::HEXAHEDRON].facets;

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

    double frame_jacobian(const mat3x3 &F) {
        // scale
        const vec3& l0 = F.col(0);
        const vec3& l1 = F.col(1);
        const vec3& l2 = F.col(2);

        // Check L_min² <= DBL_MIN => q = DBL_MAX, with DBL_MIN = 1e-30 and DBL_MAX = 1e+30
        const double l_min = std::min(std::min(l0.norm2(), l1.norm2()), l2.norm2());

        if (l_min <= 1e-30)
            return 1e30;

        double d = std::sqrt(l0.norm2() * l1.norm2() * l2.norm2());

        double scaled_jacobian = F.det() / d;

        if (scaled_jacobian > 0)
            return std::min(scaled_jacobian, 1e30);
        else 
            return std::max(scaled_jacobian, -(1e30));
    }


    // constexpr vec3 facets_of(const int v) {
    //     auto &facets = reference_cells[1].facets;
    //     auto &c2f = reference_cells[1].c2f;

    //     int cnt = 0;
    //     int idx[3]{0,0,0};
    //     for (int i = 0; i < std::size(facets); i++) {
    //         if (facets[i] == v) {
    //             idx[cnt] = i;
    //             cnt++;
    //             if (cnt == 3)
    //                 break;
    //         }
    //     }

    //     return {c2f[idx[0]], c2f[idx[1]], c2f[idx[2]]};
    // }

    // constexpr int * neighbors_verts(const int v) {

    //     auto &facets = reference_cells[1].facets;
    //     auto &c2f = reference_cells[1].c2f;

    //     int cnt = 0;
    //     int idx[3]{0,0,0};
    //     for (int i = 0; i < std::size(facets); i++) {
    //         if (facets[i] == v) {
    //             idx[cnt] = i;
    //             cnt++;
    //             if (cnt == 3)
    //                 break;
    //         }
    //     }

    //     return new int[]{
    //         reference_cells[1].vert(c2f[idx[0]], ((idx[0] % 4) + 1) % 4), 
    //         reference_cells[1].vert(c2f[idx[1]], ((idx[1] % 4) + 1) % 4), 
    //         reference_cells[1].vert(c2f[idx[2]], ((idx[2] % 4) + 1) % 4)
    //     };
    // }

    // double Hexahedron3::jacobian(int c) const {

    //     int ct = Volume::CELL_TYPE::HEXAHEDRON;
    //     // Get neightbors vertices of vertex corner
    //     const auto &neightbors = reference_cells[ct].neighbors_verts(c);
        
    //     const vec3 &pc = v[c];
    //     const vec3 &p0 = v[neightbors[0]];
    //     const vec3 &p1 = v[neightbors[1]];
    //     const vec3 &p2 = v[neightbors[2]];

    //     std::cout << "c: " << c << ", neightbors[0]:" << neightbors[0] << ",neightbors[1]:" << neightbors[1] << ",neightbors[2]:" << neightbors[2] << std::endl;
    //     std::cout << "pc:  " << pc << ", p0:" << p0 << ",p1:" << p1 << ",p2:" << p2 << std::endl;


    //     const mat3x3 FT{{p0 - pc, p1 - pc, p2 - pc}};
    //     std::cout << "--- FT ---" << std::endl << FT << std::endl;
    //     const mat3x3 F = FT.transpose();
    //     std::cout << "--- F ---" << std::endl << F << std::endl;

    //     // test compare
    //     mat<3,8> A = {{
    //         {v[0].x, v[1].x, v[3].x, v[2].x, v[4].x, v[5].x, v[7].x, v[6].x},
    //         {v[0].y, v[1].y, v[3].y, v[2].y, v[4].y, v[5].y, v[7].y, v[6].y},
    //         {v[0].z, v[1].z, v[3].z, v[2].z, v[4].z, v[5].z, v[7].z, v[6].z}
    //     }};
    //     const mat<8,3> &B = QH[c];
    //     mat3x3 J = A * B;
    //     std::cout << "--- J ---" << std::endl << J << std::endl;

    //     std::cout << "HEX VERT: " << std::endl;
    //     for (int i = 0; i < 8; i++) {
    //         std::cout << v[i] << std::endl;
    //     }

    //     return frame_jacobian(F);
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