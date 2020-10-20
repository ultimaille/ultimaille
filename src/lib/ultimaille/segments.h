#ifndef __SEGMENTS_H__
#define __SEGMENTS_H__

#include <vector>
#include <memory>
#include "geometry.h"
#include "pointset.h"

struct GenericAttributeContainer;

struct Segments { // polygonal mesh interface
    PointSet points{};
    std::vector<int> segments{};
    std::vector<std::weak_ptr<GenericAttributeContainer> > attr{};

    Segments() = default;

    int nverts() const;
    int nsegments() const;

    int create_segments(const int n);
    void resize_attrs();

    int  vert(const int s, const int lv) const;
    int &vert(const int s, const int lv)      ;
};

#endif //__SEGMENTS_H__

