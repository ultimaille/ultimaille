#include "pointset.h"
#include "attributes.h"

namespace UM {
    void PointSet::resize(const int n) {
        data->resize(n);
        resize_attrs();
    }

    int PointSet::create_points(const int n) {
        assert(n>=0);
        data->resize(size()+n);
        resize_attrs();
        return size()-n;
    }

    int PointSet::push_back(const vec3 &p) {
        data->push_back(p);
        resize_attrs();
        return size()-1;
    }

    void PointSet::resize_attrs() {
        for (auto &wp : *attr)  if (auto spt = wp.lock())
            spt->resize(size());
    }
}

