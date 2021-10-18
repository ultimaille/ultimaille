#include <algorithm>
#include "pointset.h"
#include "attributes.h"
#include "syntactic-sugar/assert.h"
#include "algebra/eigen.h"

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

    void PointSet::Util::principal_axes(vec3 &center, vec3 axes[3], double eigen_values[3]) const {
        axes[0] = {1,0,0}; // If the system is under-determined, return the trivial basis
        axes[1] = {0,1,0};
        axes[2] = {0,0,1};
        eigen_values[2] = eigen_values[1] = eigen_values[0] = 1;
        center = barycenter();
        if (ps.size()<4) return;

        mat3x3 M = {}; // covariance matrix
        for (const vec3 &p : ps)
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    M[i][j] += (p-center)[i]*(p-center)[j];
        M /= static_cast<double>(ps.size());

        vec3 eval;
        mat3x3 evec;
        eigendecompose_symmetric(M, eval, evec);

        for (int i=0; i<3; i++) {
            eigen_values[i] = eval[i];
            axes[i] = evec.col(i);
        }
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
        for (auto &wp : attr)  if (auto spt = wp.lock())
            spt->compress(old2new);
    }
}

