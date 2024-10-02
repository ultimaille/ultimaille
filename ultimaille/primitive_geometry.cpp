#include "primitive_geometry.h"
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"
#include "algebra/mat.h"
#include "volume_reference.h"
#include "volume.h"

namespace UM {
    double Segment2::distance(const vec2 &p) const {
        const double l2 = length2();
        if (l2 < 1e-15) return (p-a).norm(); // degenerate segment
        /**
         Consider the line extending the segment, parameterized as a + t (b - a).
        We find projection of point p onto the line.
        It falls where t = [(p-a) . (b-a)] / |b-a|^2
        We clamp t from [0,1] to handle points outside the segment ab.
        **/
        const double t = std::max(0., std::min(1., (p - a)*(b - a) / l2));
        const vec2 projection = a + t * (b - a);  // projection falls on the segment
        return (p-projection).norm();
    }

    vec3 Triangle2::bary_coords(vec2 G) const {
        vec3 result;
        double sum = 0;

        for (int i = 0; i < 3; i++) {
            const vec2 &A = v[(i + 1) % 3];
            const vec2 &B = v[(i + 2) % 3];

            result[i] = Triangle2(A, B, G).signed_area();
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

            result[i] = Triangle3(A, B, G).cross_product() * n;
            sum += result[i];
        }

        return result / sum;
    }

    // vec4 Tetrahedron::bary_coords(vec3 G) const {
    //     vec4 result;
    //     double sum = 0;
        
    //     double vol = volume();

    //     for (int coord = 0; coord < 4; coord++) {
    //         const vec3 &A = v[(coord + 1) % 4];
    //         const vec3 &B = v[(coord + 2) % 4];
    //         const vec3 &C = v[(coord + 3) % 4];

    //         result[coord] = std::abs(Tetrahedron(A, B, C, G).volume() / vol);
    //         sum += result[coord];
    //     }

    //     return result / sum;
    // }

    vec4 Tetrahedron::bary_coords(vec3 G) const {
    	vec4 result;
    	double sum = 0;
        
    	for(int coord=0;coord< 4;coord++) {
    		vec3 P = v[coord];
    		vec3 A = v[(coord + 1) % 4];
    		vec3 B = v[(coord + 2) % 4];
    		vec3 C = v[(coord + 3) % 4];
    		result[coord] = Tetrahedron(A, B, C, G).volume() / Tetrahedron(A, B, C, P).volume();
    		sum += result[coord];
    	}

    	result /= sum;
    	return result;
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

        return Triangle2(z0, z1, z2);
    }

    Triangle2 Triangle2::dilate(double scale) const {
        vec2 G = bary_verts();
        return Triangle2(G + scale * (v[0] - G), G + scale * (v[1] - G), G + scale * (v[2] - G));
    }

    Triangle3 Triangle3::dilate(double scale) const {
        vec3 G = bary_verts();
        return Triangle3(G + scale * (v[0] - G), G + scale * (v[1] - G), G + scale * (v[2] - G));
    }
    
    mat3x3 Triangle3::tangent_basis() const {
        mat3x3 res{v[1] - v[0], v[2] - v[0]};
        
        for (int i = 0; i < 2; i++)
            res[i].normalize();
        
        res[2] = cross(res[0], res[1]);
        res[2].normalize();
        res[1] = cross(res[2], res[0]);

        return res;
    }

    mat3x3 Triangle3::tangent_basis(vec3 first_axis) const {
        mat3x3 res{v[1] - v[0], v[2] - v[0]};
        
        for (int i = 0; i < 2; i++)
            res[i].normalize();
        
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

    mat<3,4> Tetrahedron::grad_operator() const {
        mat<4,3>  accum ;

        for (int i = 0; i < 4; i++) {
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

    vec3 Tetrahedron::grad(vec4 u) const {
        return grad_operator() * u;
    }

    double Quad3::unsigned_area() const {
        return UM::geo::unsigned_area(v, 4);
    }

    vec3 Quad3::normal() const {
        return UM::geo::normal(v, 4);
    }

    double Quad2::scaled_jacobian() const {
        // https://coreform.com/papers/verdict_quality_library.pdf
        // Note: values of DBL_MIN / DBL_MAX were found in verdict library code

        // quadratures for every quad corner
        constexpr int cverts[4][2][2] { {{0,1},{0,3}}, {{0,1},{1,2}}, {{3,2},{1,2}}, {{3,2},{0,3}} };
        
        double min_sj = std::numeric_limits<double>::max();
        for (int c = 0; c < 4; c++) { // 4 corners of the quad
            vec2 n1 = v[cverts[c][0][1]] - v[cverts[c][0][0]];
            vec2 n2 = v[cverts[c][1][1]] - v[cverts[c][1][0]];

            // Check L < DBL_MIN => q = 0 
            
            if (n1.norm() < 1.0e-30 || n2.norm() < 1.0e-30)
                return 0;

            n1.normalize();
            n2.normalize();

            min_sj = std::min(min_sj, n1.x*n2.y - n2.x*n1.y);
        }

        return min_sj;
    }

    vec3 Poly3::bary_verts() const {
        return UM::geo::bary_verts(v.data(), static_cast<int>(v.size()));
    }

    vec3 Poly3::normal() const {
        return UM::geo::normal(v.data(), static_cast<int>(v.size()));
    }

    double Poly3::unsigned_area() const {
        return UM::geo::unsigned_area(v.data(), static_cast<int>(v.size()));
    }

    vec3 Tetrahedron::bary_verts() const {
        return UM::geo::bary_verts(v, 4);
    }

    vec3 Hexahedron::bary_verts() const {
        return UM::geo::bary_verts(v, 8);
    }

    double Hexahedron::volume() const {
        const vec3 bary = bary_verts();
        double vol = 0;

        // Indexes of hex faces
        // Vertices numbering convention is very IMPORTANT (and respected in volume_reference.h)
        // as well as counter-clock wise convention for facet orientation
        const auto &indexes = reference_cells[Volume::CELL_TYPE::HEXAHEDRON].facets;

        for (int f=0; f<6; f++) {
            for (int i=0; i<4; i++) {
                const int offset = 4*f;
                vol += geo::tet_volume(
                        bary,
                        v[indexes[i + offset]],
                        v[indexes[(i+1)%4 + offset]],
                        v[indexes[(i+2)%4 + offset]]
                        )*.5;
            }
        }
        return vol;
    }

    double Hexahedron::scaled_jacobian() const {
        // TODO add extreme values condition !
        // https://coreform.com/papers/verdict_quality_library.pdf
        // Note: values of DBL_MIN / DBL_MAX were found in verdict library code

        constexpr int cverts[8][4] = { {0,1,2,4}, {1,3,0,5}, {2,0,3,6}, {3,2,1,7}, {4,6,5,0}, {5,4,7,1}, {6,7,4,2}, {7,5,6,3} };
        double min_sj = std::numeric_limits<double>::max();
        double l_min2 = std::numeric_limits<double>::max();
        for (int c = 0; c < 8; c++) { // eight corners of the cube
            vec3 n1 = v[cverts[c][1]] - v[cverts[c][0]];
            vec3 n2 = v[cverts[c][2]] - v[cverts[c][0]];
            vec3 n3 = v[cverts[c][3]] - v[cverts[c][0]];

            l_min2 = std::min(std::min(std::min(l_min2, n1.norm2()), n2.norm2()), n3.norm2());
            
            n1.normalize();
            n3.normalize();
            n2.normalize();

            min_sj = std::min(min_sj, n3 * cross(n1, n2));
        }
        { // principal axes
            vec3 n1 = v[1]-v[0] + v[3]-v[2] + v[5]-v[4] + v[7]-v[6];
            vec3 n2 = v[2]-v[0] + v[3]-v[1] + v[6]-v[4] + v[7]-v[5];
            vec3 n3 = v[4]-v[0] + v[5]-v[1] + v[6]-v[2] + v[7]-v[3];

            l_min2 = std::min(std::min(std::min(l_min2, n1.norm2()), n2.norm2()), n3.norm2());

            n1.normalize();
            n3.normalize();
            n2.normalize();

            min_sj = std::min(min_sj, n3 * cross(n1, n2));
        }

        // L_min² <= DBL_MIN => q = DBL_MAX
        if (l_min2 <= 1.0e-30)
            return 1.0e+30;

        return min_sj;
    }

    vec3 Pyramid::bary_verts() const {
        return UM::geo::bary_verts(v, 5);
    }

    double Pyramid::volume() const {
        
        const vec3 apx = apex();
        double vol = 0;

        // Compute base to apex volume
        for (int i=0; i<4; i++) {
            vol += geo::tet_volume(
                    apx,
                    v[(i+2)%4],
                    v[(i+1)%4],
                    v[i]
                    )*.5;
        }

        return vol;
    }

    vec3 Wedge::bary_verts() const {
        return UM::geo::bary_verts(v, 6);
    }

    double Wedge::volume() const {
        const vec3 bary = bary_verts();
        double vol = 0;

        // Indexes of hex faces
        // Vertices numbering convention is very IMPORTANT (and respected in volume_reference.h)
        // as well as counter-clock wise convention for facet orientation
        auto rc = reference_cells[Volume::CELL_TYPE::WEDGE];
        const auto &indexes = rc.facets;

        vol += geo::tet_volume(
            bary,
            v[indexes[0]],
            v[indexes[1]],
            v[indexes[2]]
        );

        vol += geo::tet_volume(
            bary,
            v[indexes[3]],
            v[indexes[4]],
            v[indexes[5]]
        );

        for (int f=0; f<rc.nf-2; f++) {

            for (int i=0; i<4; i++) {
                const int offset = 6+4*f;
                vol += geo::tet_volume(
                        bary,
                        v[indexes[i + offset]],
                        v[indexes[(i+1)%4 + offset]],
                        v[indexes[(i+2)%4 + offset]]
                        )*.5;
            }
        }
        return vol;
    }

}
