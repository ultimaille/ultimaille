#include <vector>
#include <memory>
#include "surface.h"
#include "attributes.h"

template struct FacetAttribute<int>;
template struct FacetAttribute<double>;

template <typename T> FacetAttribute<T>::FacetAttribute(Surface &m) : ac(new AttributeContainer<T>(m.nfacets())) {
    m.attr_facets.push_back(ac);
}

