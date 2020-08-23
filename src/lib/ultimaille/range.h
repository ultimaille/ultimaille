#ifndef __RANGE_H__
#define __RANGE_H__

#include "surface.h"

struct range {
    const int end_;

    range(const int end) : end_(end) {}

    struct iterator {
        iterator &operator++() { ++value_; return *this; }
        const int &operator*() const { return value_; }
        bool operator!=(const iterator& rhs) const { return value_ != rhs.value_; }
        int value_{};
    };

    iterator begin() const { return iterator{0};    }
    iterator end()   const { return iterator{end_}; }
};

inline range facets(Surface &m) {
    return {m.nfacets()};
}

inline range corners(Surface &m) {
    return {m.ncorners()};
}

struct facet_vertices {
    Surface &m_;
    const int facet_;
    facet_vertices(Surface &m, const int facet) : m_(m), facet_(facet) {}

    struct iterator {
        iterator &operator++() { ++value_; return *this; }
        int &operator*() const { return m_.vert(facet_, value_); }
        bool operator!=(const iterator& rhs) const { return value_ != rhs.value_; }
        Surface &m_;
        const int facet_;
        int value_;
    };

    iterator begin() const { return iterator{m_, facet_, 0};    }
    iterator end()   const { return iterator{m_, facet_, m_.facet_size(facet_)}; }
};


#endif // __RANGE_H__

