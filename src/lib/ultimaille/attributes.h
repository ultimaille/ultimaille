#ifndef __ATTRIBUTES_H__
#define __ATTRIBUTES_H__
#include <vector>
#include <memory>
#include "surface.h"

struct GenericAttributeContainer {
    virtual void resize(const int n) = 0;
    virtual ~GenericAttributeContainer() = default;
};

template <typename T> struct AttributeContainer : GenericAttributeContainer {
    AttributeContainer(const int n) : data(n) {}
    void resize(const int n) { data.resize(n); }
    std::vector<T> data;
};

template <typename T> struct FacetAttribute {
    FacetAttribute(Surface &m) : ptr(new AttributeContainer<T>(m.nfacets())) {
        m.attr_facets.push_back(ptr);
    }
          T& operator[](const int i)       { return std::dynamic_pointer_cast<AttributeContainer<T> >(ptr)->data[i]; }
    const T& operator[](const int i) const { return std::dynamic_pointer_cast<AttributeContainer<T> >(ptr)->data[i]; }

    std::shared_ptr<GenericAttributeContainer> ptr;
};

#endif //__ATTRIBUTES_H__

