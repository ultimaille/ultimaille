#ifndef __POLYLINE_H__
#define __POLYLINE_H__

#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "pointset.h"
#include "syntactic-sugar/assert.h"

namespace UM {
    struct GenericAttributeContainer;

    struct PolyLine {
        PointSet points{};
        std::vector<int> segments{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr{};

        int nverts() const;
        int nsegments() const;

        void compress_attrs(const std::vector<bool> &segments_to_kill);
        void delete_vertices(const std::vector<bool> &to_kill);
        void delete_segments(const std::vector<bool> &to_kill);
        int create_segments(const int n);
        void resize_attrs();

        int  vert(const int s, const int lv) const;
        int &vert(const int s, const int lv)      ;

        virtual void clear() {
            points   = {};
            attr     = {};
            segments = {};
        }

        PolyLine() : util(*this) {}
        PolyLine(const PolyLine& m) : util(*this) {
            um_assert(!m.points.size() && !m.segments.size());
        }
        PolyLine& operator=(const PolyLine& m) {
            clear();
            um_assert(!m.points.size() && !m.segments.size());
            return *this;
        }

        struct Util {
            Util(const PolyLine &mesh) : m(mesh) {}
            const PolyLine &m;
        } util;
    };
}

#endif //__POLYLINE_H__

