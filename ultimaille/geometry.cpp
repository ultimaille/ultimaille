#include "geometry.h"
#include <initializer_list>
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
        double a = 0;
        vec3 G = bary_verts();
        for (int lv=0; lv<4; lv++) {
            const vec3 &A = v[lv];
            const vec3 &B = v[(lv+1)%4];
            a += UM::unsigned_area(G, A, B);
        }
        return a;
    }

}