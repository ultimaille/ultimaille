#include "geometry.h"
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"
#include "algebra/mat.h"
#include "volume_reference.h"

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
        
        for (int i = 0; i < 4; i++) {
            const vec3 &P = v[i];
            const vec3 &A = v[(i + 1) % 4];
            const vec3 &B = v[(i + 2) % 4];
            const vec3 &C = v[(i + 3) % 4];

            result[i] = Tetrahedron3{{A, B, C, G}}.volume() / Tetrahedron3{{A, B, C, P}}.volume();
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
        
        for (int d = 0; d < 3; d++) {
            vec3 p = v[d];
            vec3 a = v[(d + 1) % 3];
            vec3 b = v[(d + 2) % 3];
            vec3 n = cross(b - a, tr_normal);
            n.normalize();
            double scale = (n * (p - a));
            accum[d] = n / scale;
        }

        return accum.transpose();
    }

    vec3 Triangle3::grad(vec3 u) const {
        return grad_operator() * u;
    }

	mat<2,3> Triangle2::grad_operator() const {
		mat<3,2>  accum;

        for (int d = 0; d < 3; d++) {
			vec2 P = v[d];
			vec2 A = v[(d + 1) % 3];
			vec2 B = v[(d + 2) % 3];

			vec2 n((B - A)[1], -(B - A)[0]);
			n.normalize();
			double scale = (n * (P - A));
			accum[d] = n / scale;
		}
		return accum.transpose();
	}

	vec2 Triangle2::grad(vec3 u) const {
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

    double Hexahedron3::jacobian(int c) const {

        // Convert verdict volume reference to ultimaille volume reference (https://coreform.com/papers/verdict_quality_library.pdf, p.79)
        // Note: Only need to invert points at index 2,3 and 6,7
        mat<3,8> A = {{
            {v[0].x, v[1].x, v[3].x, v[2].x, v[4].x, v[5].x, v[7].x, v[6].x},
            {v[0].y, v[1].y, v[3].y, v[2].y, v[4].y, v[5].y, v[7].y, v[6].y},
            {v[0].z, v[1].z, v[3].z, v[2].z, v[4].z, v[5].z, v[7].z, v[6].z}
        }};

        for (int lc = 0; lc < 8; lc++)
            for (int f = 0; f < 3; f++) {
                const int he = reference_cells[1].corner(f, lc);
                // std::cout << lc << "," << f << " -> " << reference_cells[1].from(he) << std::endl;
            }

        // p1-p0, p2-p0, p4-p0
        // 


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