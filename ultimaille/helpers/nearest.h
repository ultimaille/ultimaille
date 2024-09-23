#ifndef __NEAREST_H__
#define __NEAREST_H__

#include <vector>
#include "ultimaille/algebra/vec.h"
#include "ultimaille/helpers/hboxes.h"
#include "ultimaille/surface.h"

namespace UM {

    struct PointOnMesh {
        inline operator vec3() const { return p; }
        Surface::Facet f;
        vec3 p;
    };

    struct NearestPointOnMesh : HBoxes<3> {
        NearestPointOnMesh(const Triangles &m);
        PointOnMesh query(vec3 p);

        Triangles &m; // TODO convert it to const ref
        std::vector<BBox3> bboxes;
    };

}

#endif //__NEAREST_H__

