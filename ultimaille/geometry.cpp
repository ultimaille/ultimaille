#include "geometry.h"
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"


namespace UM {

    void Triangle3::project(vec2& z0, vec2& z1, vec2& z2) const {
        const vec3 &A = v[0];
        const vec3 &B = v[1];
        const vec3 &C = v[2];

        vec3 X = (B - A).normalized(); // construct an orthonormal 3d basis
        vec3 Z = cross(X, C - A).normalized();
        vec3 Y = cross(Z, X);

        z0 = vec2(0,0); // project the triangle to the 2d basis (X,Y)
        z1 = vec2((B - A).norm(), 0);
        z2 = vec2((C - A)*X, (C - A)*Y);
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

    vec3 Poly3::bary_verts() const {
        return UM::bary_verts(v.data(), static_cast<const int>(v.size()));
    }

    vec3 Tetrahedron3::bary_verts() const {
        return UM::bary_verts(v, 4);
    }

    vec3 Hexahedron3::bary_verts() const {
        return UM::bary_verts(v, 8);
    }

    vec3 Pyramid3::bary_verts() const {
        return UM::bary_verts(v, 5);
    }

    double Hexahedron3::volume() const {
        const vec3 bary = bary_verts();
        double vol = 0;

        // Indexes of hex faces
        // Vertices numbering convention (from geogram) is very IMPORTANT (and respected)
        // as well as counter-clock wise convention for facet orientation
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

    double Pyramid3::volume() const {
        // TODO fill
        return 0;
    }

    vec3 Poly3::normal() const {
        return UM::normal(v.data(), static_cast<const int>(v.size()));
    }

    double Poly3::unsigned_area() const {
        return UM::unsigned_area(v.data(), static_cast<const int>(v.size()));
    }

}