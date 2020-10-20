#ifndef __POLYLINE_H__
#define __POLYLINE_H__

#include <vector>
#include <memory>
#include "geometry.h"
#include "pointset.h"

struct GenericAttributeContainer;

struct PolyLine { // polygonal mesh interface
    PointSet points{};
    std::vector<int> segments{};
    std::vector<std::weak_ptr<GenericAttributeContainer> > attr{};

    PolyLine() = default;

    int nverts() const;
    int nsegments() const;

    int create_segments(const int n);
    void resize_attrs();

    int  vert(const int s, const int lv) const;
    int &vert(const int s, const int lv)      ;
};

#endif //__POLYLINE_H__

