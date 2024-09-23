#include <limits>
#include <queue>
#include "ultimaille/primitive_geometry.h"
#include "ultimaille/helpers/nearest.h"

namespace UM {

    BVHTriangles::BVHTriangles(const Triangles &mesh) : m(const_cast<Triangles &>(mesh)) {
        std::vector<BBox3> bboxes(m.nfacets());
        for (int f=0; f<m.nfacets(); f++)               // create boxes bounding
            for (int lv=0; lv<m.facet_size(f); lv++)    // individual faces
                bboxes[f].add(m.points[m.vert(f, lv)]);
        init(bboxes);                                   // create the bounding volume hierarchy
    }

    inline double dist_segment(double a, double b, double x) {
            return x < a ? a-x : (x > b ? x-b : 0.);
    }

    inline double dist2_box(const BBox3 &box, const vec3 &p) {
        return vec3(
                dist_segment(box.min.x, box.max.x, p.x),
                dist_segment(box.min.y, box.max.y, p.y),
                dist_segment(box.min.z, box.max.z, p.z)
                ).norm2();
    }

    PointOnMesh BVHTriangles::nearest_point(vec3 p) {
        double best_dist2 = std::numeric_limits<double>::max();
        PointOnMesh best_point = {{m, 0}, {0,0,0}};

        using QEl = std::pair<double, int>;
        std::priority_queue<QEl, std::vector<QEl>, std::greater<QEl>> Q;
        Q.emplace(0., 0);

        while (!Q.empty() && Q.top().first < best_dist2) {
            const int node = Q.top().second; Q.pop();
            const int leaves = tree.size()  - m.nfacets();               // start offset for the leaves of the hierarchy
            const int beg = 2*node + 1;                                  // start offset for the children nodes
            const int end = std::min(                                    //   end offset for the children nodes
                    2*node + 3,
                    static_cast<int>(tree.size())
                    );

            for (int son = beg; son<end; son++) {                        // iterate through children boxes
                if (son < leaves)                                        // if it is not a leaf, place it in the priority queue
                    Q.emplace(dist2_box(tree[son], p), son);
                else {
                    Surface::Facet f = {m, tree_pos_to_org[son-leaves]}; // for the leaves we can directly compute
                    vec3 nearest = Triangle3(f).nearest_point(p);        // the nearest point and compare it to the current best
                    double dist2 = (p-nearest).norm2();
                    if (best_dist2 > dist2) {
                        best_dist2 = dist2;
                        best_point = {f, nearest};
                    }
                }
            }
        }
        return best_point;
    }

}

