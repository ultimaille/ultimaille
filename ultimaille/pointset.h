#ifndef __POINTSET_H__
#define __POINTSET_H__
#include <vector>
#include <tuple>
#include <memory>
#include "algebra/vec.h"
#include "algebra/mat.h"
#include "helpers/hboxes.h"

namespace UM {
    struct ContainerBase;

    struct PointSet {
        PointSet() : data(new std::vector<vec3>()) {}
        PointSet(std::shared_ptr<std::vector<vec3> > ext) : data(ext) {}

        PointSet(const PointSet &p)            = default; // We need to be able to share point sets, therefore we allow copying of the pointer
        PointSet(PointSet &&p)                 = default; // N.B. attr pointers are also copied, but this should not have 
        PointSet& operator=(const PointSet& p) = default; // serious consequences since use_count()==1 is asserted for modification

        int size() const { return data->size(); }
        vec3& operator[](const int i) { return data->at(i); }
        const vec3& operator[](const int i) const { return data->at(i); }
        int use_count() { return data.use_count(); }

        void resize(const int n);
        int push_back(const vec3 &p);

        void delete_points(const std::vector<bool> &to_kill);
        void delete_points(const std::vector<bool> &to_kill, std::vector<int> &old2new); // TODO: remove old2new
        int create_points(const int n);

        using       iterator = std::vector<vec3>::iterator;
        using const_iterator = std::vector<vec3>::const_iterator;

        iterator begin() { return data->begin(); }
        iterator end()   { return data->end();   }
        const_iterator begin() const { return data->begin(); }
        const_iterator end()   const { return data->end();   }

        void resize_attrs();
        void compress_attrs(const std::vector<int> &old2new);

        std::shared_ptr<std::vector<vec3> > data;
        std::vector<std::weak_ptr<ContainerBase> > attr = {};
    };
}

#endif //__POINTSET_H__

