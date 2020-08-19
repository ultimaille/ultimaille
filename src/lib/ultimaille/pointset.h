#ifndef __POINTSET_H__
#define __POINTSET_H__
#include <vector>
#include <memory>
#include "geometry.h"

struct GenericAttributeContainer;

struct PointSet {
    PointSet() : data(new std::vector<vec3>()), attr() {}
    PointSet(std::shared_ptr<std::vector<vec3> > ext) : data(ext), attr() {}

          vec3& operator[](const int i)       { return data->at(i); }
    const vec3& operator[](const int i) const { return data->at(i); }

    int size() const { return data->size(); }

    void bbox(vec3 &min, vec3 &max) {
        min = max = data->at(0);
        for (vec3 const &p : *data) {
            for (int j=0; j<3; j++) {
                min[j] = std::min(min[j], p[j]);
                max[j] = std::max(max[j], p[j]);
            }
        }
    }

    void push_back(const vec3 &p) {
        assert(1==data.use_count());
        data->push_back(p);
        resize_attrs();
    }

    void delete_points(const std::vector<bool> &to_kill, std::vector<int> &old2new) {
        assert(1==data.use_count());
        assert(to_kill.size()==(size_t)size());
        old2new = std::vector<int>(size(),  -1);

        int new_nb_pts = 0;
        for (int v=0; v<size(); v++) {
            if (to_kill[v]) continue;
            data->at(new_nb_pts) = data->at(v);
            old2new[v] = new_nb_pts++;
        }
        data->resize(new_nb_pts);
        compress_attrs(old2new);
    }

    using       iterator = std::vector<vec3>::iterator;
    using const_iterator = std::vector<vec3>::const_iterator;

          iterator begin()       { return data->begin(); }
    const_iterator begin() const { return data->begin(); }
          iterator end()       { return data->end(); }
    const_iterator end() const { return data->end(); }

    void resize_attrs();
    void compress_attrs(const std::vector<int> &old2new);

    std::shared_ptr<std::vector<vec3> > data;
    std::vector<std::weak_ptr<GenericAttributeContainer> > attr;
};

#endif //__POINTSET_H__

