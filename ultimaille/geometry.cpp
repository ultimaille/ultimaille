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

    vec3 Poly3::normal() const {
        return UM::normal(v.data(), static_cast<const int>(v.size()));
    }

    double Poly3::unsigned_area() const {
        return UM::unsigned_area(v.data(), static_cast<const int>(v.size()));
    }

}