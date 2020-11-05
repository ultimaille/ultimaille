#ifndef __HBOXES_H__
#define __HBOXES_H__

#include <vector>
#include <limits>
#include "geometry.h"

namespace UM {
    struct BBox3 {
        BBox3() {
            min = { std::numeric_limits<double>::max(),  std::numeric_limits<double>::max(),  std::numeric_limits<double>::max()};
            max = {-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
        }

        void add(const BBox3 &b) {
            if (b.empty()) return;
            add(b.min);
            add(b.max);
        }

        void add(const vec3 &p) {
            for (int d=0; d<3; d++) {
                min[d] = std::min<double>(min[d], p[d]);
                max[d] = std::max<double>(max[d], p[d]);
            }
        }

        bool empty() const {
            for (int d=0; d<3; d++)
                if (max[d] < min[d])
                    return true;
            return false;
        }

        bool intersect(const BBox3 &b) const {
            for (int d=0; d<3; d++)
                if (min[d] > b.max[d] || max[d] < b.min[d])
                    return false;
            return true;
        }

        bool contains(const vec3 &v) const {
            for (int d=0; d<3; d++)
                if (min[d] > v[d] || max[d] < v[d])
                    return false;
            return true;
        }

        void dilate(double eps) {
            min = min - vec3(eps, eps, eps);
            max = max + vec3(eps, eps, eps);
        }

        vec3 center() const {
            return (min + max)*.5;
        }

        vec3 min;
        vec3 max;
    };

    struct HBoxes {
        HBoxes(std::vector<BBox3> const &inboxes);

        void sort(std::vector<vec3> &G, int org, int dest) const;
        void intersect(BBox3 const &b, std::vector<int> &primitives, int node = 0) const;

        int offset;
        mutable std::vector<int> tree_pos_to_org;
        std::vector<BBox3> tree;
    };
}

#endif //__HBOXES_H__

