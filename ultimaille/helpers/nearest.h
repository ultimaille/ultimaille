#ifndef __NEAREST_H__
#define __NEAREST_H__

#include "ultimaille/algebra/vec.h"
#include "ultimaille/helpers/hboxes.h"
#include "ultimaille/surface.h"

namespace UM {

    struct PointOnMesh {
        inline operator vec3() const { return p; }
        Surface::Facet f;
        vec3 p;
    };

    struct BVHTriangles : HBoxes<3> {
        BVHTriangles(const Triangles &m);
        PointOnMesh nearest_point(vec3 p);

        Triangles &m; // TODO convert it to const ref
    };

}

#endif //__NEAREST_H__

