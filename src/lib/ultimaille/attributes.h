#ifndef __ATTRIBUTES_H__
#define __ATTRIBUTES_H__
#include <vector>
#include <memory>
#include "surface.h"
#include "permutation.h"

struct GenericAttributeContainer {
    virtual void resize(const int n) = 0;
    virtual void permute(Permutation &perm) = 0;
    virtual ~GenericAttributeContainer() = default;
};

template <typename T> struct AttributeContainer : GenericAttributeContainer {
    AttributeContainer(const int n) : data(n) {}
    void resize(const int n) { data.resize(n); }
    void permute(Permutation &perm) {
        perm.apply(data);
    }
    std::vector<T> data;
};

template <typename T> struct GenericAttribute {
    GenericAttribute(int size) : ptr(new AttributeContainer<T>(size)) {}
    T& operator[](const int i)       { return std::dynamic_pointer_cast<AttributeContainer<T> >(ptr)->data[i]; }
    T  operator[](const int i) const { return std::dynamic_pointer_cast<AttributeContainer<T> >(ptr)->data[i]; }
    std::shared_ptr<GenericAttributeContainer> ptr;
};

template <typename T> struct FacetAttribute : GenericAttribute<T> {
    FacetAttribute(Surface &m) : GenericAttribute<T>(m.nfacets())  {
        m.attr_facets.push_back(this->ptr);
    }
};

#endif //__ATTRIBUTES_H__

