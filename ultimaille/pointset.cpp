#include <algorithm>
#include "pointset.h"
#include "attributes.h"
#include "syntactic-sugar/assert.h"
#include "algebra/eigen.h"
#include "algebra/covariance.h"

namespace UM {
    BBox3 PointSet::Util::bbox() const {
        BBox3 bbox;
        for (vec3 const &p : ps)
            bbox.add(p);
        return bbox;
    }

    vec3 PointSet::Util::barycenter() const {
        vec3 ave = {0, 0, 0};
        for (const vec3 &p : ps)
            ave += p;
        return ave / static_cast<double>(ps.size());
    }

    std::tuple<mat3x3,vec3,vec3> PointSet::Util::principal_axes() const {
        PointSetCovariance cov(*ps.data);
        vec3 eval;
        mat3x3 evec;
        eigendecompose_symmetric(cov.cov, eval, evec);
        if (ps.size()<4) return {mat3x3::identity(), {1.,1.,1.}, cov.center}; // If the system is under-determined, return the trivial basis
        return { evec, eval, cov.center };
    }

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

    void PointSet::delete_points(const std::vector<bool> &to_kill, std::vector<int> &old2new) {
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

    void PointSet::resize_attrs() {
        um_assert(1==data.use_count());
        for (auto &wp : attr)  if (auto spt = wp.lock())
            spt->resize(size());
    }

    void PointSet::compress_attrs(const std::vector<int> &old2new) {
        um_assert(1==data.use_count());
        std::erase_if(attr, [](std::weak_ptr<GenericAttributeContainer> ptr) { return ptr.lock()==nullptr; }); // remove dead attributes
        for (auto &wp : attr) { // compress attributes
            auto spt = wp.lock();
            assert(spt!=nullptr);
            spt->compress(old2new);
        }
    }
}

