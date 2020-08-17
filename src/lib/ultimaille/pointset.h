#ifndef __POINTSET_H__
#define __POINTSET_H__
#include <vector>
#include <memory>
#include "geometry.h"

struct PointSet {
    PointSet() : data(new std::vector<vec3>()) {}
    PointSet(std::shared_ptr<std::vector<vec3> > ext) : data(ext) {}

          vec3& operator[](const int i)       { return data->at(i); }
    const vec3& operator[](const int i) const { return data->at(i); }

    int size() const { return data->size(); }

    void get_bbox(vec3 &min, vec3 &max) {
        min = max = data->at(0);
        for (vec3 const &p : *data) {
            for (int j=0; j<3; j++) {
                min[j] = std::min(min[j], p[j]);
                max[j] = std::max(max[j], p[j]);
            }
        }
    }

    void push_back(const vec3 &p) { data->push_back(p); }

    using       iterator = std::vector<vec3>::iterator;
    using const_iterator = std::vector<vec3>::const_iterator;

          iterator begin()       { return data->begin(); }
    const_iterator begin() const { return data->begin(); }
          iterator end()       { return data->end(); }
    const_iterator end() const { return data->end(); }

    std::shared_ptr<std::vector<vec3> > data;
};

#endif //__POINTSET_H__

