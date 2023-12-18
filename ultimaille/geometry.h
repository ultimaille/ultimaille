#ifndef __GEOM_H__
#define __GEOM_H__
#include <initializer_list>
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"

namespace UM {

    enum CELL_TYPE { TRI=0, QUAD=1, POLY=3 };

    struct Triangle3 {
        vec3 v[3] = {};
        inline vec3 normal() const;
        inline vec3 bary_verts() const;
    };

    struct Quad3 {
        vec3 v[4] = {};
        inline vec3 normal() const;
        inline vec3 bary_verts() const;
    };

    struct Poly3 {
        std::vector<vec3> v = {};
    };

    inline vec3 Triangle3::normal() const {
        return cross(v[1]-v[0], v[2]-v[0]).normalized();
    }

    inline vec3 Triangle3::bary_verts() const {
        return (v[0] + v[1] + v[2]) / 3;
    }

    // Copy past of Utils !
    inline vec3 Quad3::normal() const {
        vec3 res = {0, 0, 0};
        vec3 bary = bary_verts();
        for (int lv=0; lv<4; lv++)
            res += cross(
                v[lv]-bary,
                v[(lv+1)%4]-bary
            );
        return res.normalized();
    }

    // Copy past of Utils !
    inline vec3 Quad3::bary_verts() const {
        vec3 ave = {0, 0, 0};
        const int nbv = 4;
        for (int lv=0; lv<nbv; lv++)
            ave += v[lv];
        return ave / static_cast<double>(nbv);
    }

}

#endif //__GEOM_H__