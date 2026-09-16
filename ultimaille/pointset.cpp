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

    void PointSet::compress_attrs(const std::vector<int> &old2new) {
        std::erase_if(*attr, [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; }); // remove dead attributes
        for (auto &wp : *attr) { // compress attributes
            auto spt = wp.lock();
            assert(spt!=nullptr);
            spt->compress(old2new);
        }
    }
}

