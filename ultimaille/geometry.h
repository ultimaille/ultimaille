#ifndef __GEOM_H__
#define __GEOM_H__
#include <initializer_list>
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "syntactic-sugar/assert.h"

namespace UM {

    struct Triangle3d {
        vec3 v[3] = {};
        inline vec3 normal() const;
        inline vec3 bary_verts() const;
    };

    inline vec3 Triangle3d::normal() const {
        return cross(v[1]-v[0], v[2]-v[0]).normalized();
    }

    inline vec3 Triangle3d::bary_verts() const {
        return (v[0] + v[1] + v[2]) / 3;
    }

    struct Quad3d {
        vec3 v[4] = {};
        inline vec3 normal() const;
        inline vec3 bary_verts() const;
    };

    // Copy past of Utils !
    inline vec3 Quad3d::bary_verts() const {
        vec3 ave = {0, 0, 0};
        const int nbv = 4;
        for (int lv=0; lv<nbv; lv++)
            ave += v[lv];
        return ave / static_cast<double>(nbv);
    }

    // Copy past of Utils !
    inline vec3 Quad3d::normal() const {
        vec3 res = {0, 0, 0};
        vec3 bary = bary_verts();
        for (int lv=0; lv<4; lv++)
            res += cross(
                    v[lv]-bary,
                    v[(lv+1)%4]-bary
                    );
        return res.normalized();
    }

    struct Poly3d {
        std::vector<vec3> v = {};
    };
    
}

#endif //__GEOM_H__