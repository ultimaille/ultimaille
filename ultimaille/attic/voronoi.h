#ifndef __VORONOI_H__
#define __VORONOI_H__

#include <vector>
#include "knn.h"

namespace UM {
    struct ConvexCell {
        enum Status {
            vertex_overflow = 1,
            inconsistent_boundary = 2,
            security_radius_not_reached = 3,
            success = 4,
            empty_cell = 6
        };

        ConvexCell(const vec2 min, const vec2 max);

        void clip_by_line(const vec3 eqn);
        void export_verts(std::vector<vec2> &verts) const;
        vec3 vertex(const int t) const;

        // N.B. clip[], tr[] and head suffice to compute all; vpos[], nb_v and nb_t are here for performance reasons

        static constexpr int _MAX_P_ = 256;
        vec3 clip[_MAX_P_];  // clipping half-planes
        int    tr[_MAX_P_];  // memory pool for chained lists of triangles
        vec3 vpos[_MAX_P_];  // precomputed convex cell vertices (corresponds to tr[])

        int head;
        int nb_v; // number of clipping lines
        int nb_t; // number of convex cell vertices

        Status status;
    };

    bool voronoi_cell(const KNN<2> &knn, ConvexCell &cc, std::vector<vec2> &verts,  const std::vector<vec2> &pts, const int seed, const int k);
    //void voronoi_diagram(const std::vector<vec2> &pts, Mesh &voronoi, double colocate_tolerance);
}

#endif // __VORONOI_H__

