#ifndef __GEOM_H__
#define __GEOM_H__
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"

namespace UM {

    inline double unsigned_area(const vec3 &A, const vec3 &B, const vec3 &C) {
        return 0.5*cross(B-A, C-A).norm();
    }

    inline vec3 bary_verts(const vec3 v[], const int nbv) {
        vec3 ave = {0, 0, 0};
        for (int lv=0; lv<nbv; lv++)
            ave += v[lv];
        return ave / static_cast<double>(nbv);
    }

    inline vec3 normal(const vec3 v[], const int nbv) {
        vec3 res = {0, 0, 0};
        vec3 bary = bary_verts(v, nbv);
        
        for (int lv=0; lv<nbv; lv++)
            res += cross(
                v[lv]-bary,
                v[(lv+1)%nbv]-bary
            );
        return res.normalized();
    }

    inline double unsigned_area(const vec3 v[], const int nbv) {
        double a = 0;
        vec3 G = UM::bary_verts(v, nbv);
        for (int lv=0; lv<nbv; lv++) {
            const vec3 &A = v[lv];
            const vec3 &B = v[(lv+1)%nbv];
            a += UM::unsigned_area(G, A, B);
        }
        return a;
    }

    struct Triangle2 {
        vec2 v[3] = {};
        inline vec2 bary_verts() const;
        inline double signed_area() const;
    };

    inline vec2 Triangle2::bary_verts() const {
        return (v[0] + v[1] + v[2]) / 3;
    }

    inline double Triangle2::signed_area() const {
        const vec2 &A = v[0];
        const vec2 &B = v[1];
        const vec2 &C = v[2];
        return .5*((B.y-A.y)*(B.x+A.x) + (C.y-B.y)*(C.x+B.x) + (A.y-C.y)*(A.x+C.x));
    }

    struct Triangle3 {
        vec3 v[3] = {};
        inline vec3 normal() const;
        inline vec3 bary_verts() const;
        inline double unsigned_area() const;
        Triangle2 project() const;
    };

    inline vec3 Triangle3::normal() const {
        return cross(v[1]-v[0], v[2]-v[0]).normalized();
    }

    inline vec3 Triangle3::bary_verts() const {
        return (v[0] + v[1] + v[2]) / 3;
    }

    inline double Triangle3::unsigned_area() const {
        return UM::unsigned_area(v[0], v[1], v[2]);
    }

    struct Quad3 {
        vec3 v[4] = {};
        vec3 normal() const;
        vec3 bary_verts() const;
        double unsigned_area() const;
    };

    struct Poly3 {
        std::vector<vec3> v = {};
        vec3 normal() const;
        vec3 bary_verts() const;
        double unsigned_area() const;
    };

    struct Tetrahedron3 {
        vec3 v[4] = {};
        vec3 bary_verts() const;
        inline double volume() const;
    };

    // signed volume for a tet with vertices (A,B,C,D) such that (AB, AC, AD) form a right hand basis
    inline double tet_volume(const vec3 &A, const vec3 &B, const vec3 &C, const vec3 &D) {
        return ((B-A)*cross(C-A,D-A))/6.;
    }

    inline double Tetrahedron3::volume() const {
        return tet_volume(
                v[0],
                v[1],
                v[2],
                v[3]
            );
    }

    struct Hexahedron3 {
        vec3 v[8] = {};
        double volume() const;
        vec3 bary_verts() const;
        // TODO lbinria: add ?
        double jacobian(int c) const;
        double jacobian2(int c) const;
        double scaled_jacobian() const;
    };

    struct Pyramid3 {
        vec3 v[5] = {};
        double volume() const;
        vec3 bary_verts() const;
        inline Quad3 base() const;
        inline const vec3 apex() const;
    };

    inline Quad3 Pyramid3::base() const {
        return Quad3{{v[0], v[1], v[2], v[3]}};
    }

    inline const vec3 Pyramid3::apex() const {
        return v[4];
    }

}

#endif //__GEOM_H__