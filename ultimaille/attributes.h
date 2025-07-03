#ifndef __ATTRIBUTES_H__
#define __ATTRIBUTES_H__
#include <vector>
#include <memory>
#include <cassert>
#include "syntactic-sugar/assert.h"
#include "pointset.h"
//#include "polyline.h"
//#include "surface.h"
//#include "volume.h"

namespace UM {
//    struct PointSet;
    struct PolyLine;
    struct Surface;
    struct Volume;

    struct SurfaceAttributes;
    struct PointSetAttributes;
    struct PolyLineAttributes;
    struct VolumeAttributes;

    struct ContainerBase {
        virtual ~ContainerBase() = default;
        virtual void resize(const int n) = 0;
        virtual void compress(const std::vector<int> &old2new) = 0;
    };

    struct AttributeBase {
        enum TYPE { GENERIC=-1, POINTS=0, EDGES=1, FACETS=2, CORNERS=3, CELLS=4, CELLFACETS=5, CELLCORNERS=6 };
        virtual ~AttributeBase() = default;
        virtual TYPE kind() const { return GENERIC; }
        virtual std::shared_ptr<ContainerBase> get_ptr() const = 0;
    };

    template <typename T> struct AttributeContainer : ContainerBase {
        AttributeContainer(const int n, const T def = T()) : data(n, def), default_value(def) {}
        void resize(const int n) { data.resize(n, default_value); }
        void compress(const std::vector<int> &old2new) { // NB: old2new is not a permutation!
            assert(old2new.size()==data.size());
            int cnt = 0;
            for (int i=0; i<static_cast<int>(old2new.size()); i++) {
                if (old2new[i]<0) continue;
                data[old2new[i]] = data[i];
                cnt++;
            }
            resize(cnt);
        }
        std::vector<T> data;
        T default_value;
    };

    template <typename T> struct GenericAttribute : AttributeBase {
//  using value_type = T;
        GenericAttribute() : ptr(nullptr) {}
        GenericAttribute(int size, const T def = T()) : ptr(new AttributeContainer<T>(size, def)) {}
        GenericAttribute(std::shared_ptr<AttributeContainer<T> > p) : ptr(p) {}
        GenericAttribute(const GenericAttribute<T>& rhs) = delete;
        GenericAttribute<T>& operator=(const GenericAttribute<T>& rhs) = delete;
              T& operator[](const int i)       { return ptr->data[i]; }
        const T& operator[](const int i) const { return ptr->data[i]; }
        void fill(T value) {
            if (ptr) std::fill(ptr->data.begin(), ptr->data.end(), value);
        }
        virtual std::shared_ptr<ContainerBase> get_ptr() const { return ptr; }
        std::shared_ptr<AttributeContainer<T> > ptr;
    };

    template <> struct GenericAttribute<bool> : AttributeBase {
//  using value_type = bool;
        GenericAttribute() : ptr(nullptr) {}
        GenericAttribute(int size, const bool def = false) : ptr(new AttributeContainer<bool>(size, def)) {}
        GenericAttribute(std::shared_ptr<AttributeContainer<bool> > p) : ptr(p) {}
        GenericAttribute(const GenericAttribute<bool>& rhs) = delete;
        GenericAttribute<bool>& operator=(const GenericAttribute<bool>& rhs) = delete;

        struct ConstBoolAttributeAccessor {
            ConstBoolAttributeAccessor(const GenericAttribute<bool>& attribute, int index) : attribute(&attribute), index(index) {}
            operator bool() const {
                return attribute->ptr->data[index];
            }
            const GenericAttribute<bool>* attribute;
            const int index;
        };

        struct BoolAttributeAccessor {
            BoolAttributeAccessor(GenericAttribute<bool>& attribute, int index) : attribute(&attribute), index(index) { }
            BoolAttributeAccessor(const BoolAttributeAccessor& rhs) : attribute(rhs.attribute), index(rhs.index) { }

            operator bool() const {
                return attribute->ptr->data[index];
            }

            BoolAttributeAccessor& operator=(bool x) {
                attribute->ptr->data[index] = x;
                return *this;
            }

            BoolAttributeAccessor& operator=(const BoolAttributeAccessor& rhs) {
                if (&rhs != this)
                    attribute->ptr->data[index] = rhs.attribute->ptr->data[rhs.index];
                return *this;
            }

            BoolAttributeAccessor& operator=(const ConstBoolAttributeAccessor& rhs) {
                attribute->ptr->data[index] = rhs.attribute->ptr->data[rhs.index];
                return *this;
            }

            GenericAttribute<bool>* attribute;
            int index;
        };

        BoolAttributeAccessor operator[](const int i) {
            return BoolAttributeAccessor(*this, i);
        }

        ConstBoolAttributeAccessor operator[](const int i) const {
            return ConstBoolAttributeAccessor(*this, i);
        }

        void fill(bool value) {
            if (ptr) std::fill(ptr->data.begin(), ptr->data.end(), value);
        }

        virtual std::shared_ptr<ContainerBase> get_ptr() const { return ptr; }

        std::shared_ptr<AttributeContainer<bool> > ptr;
    };

    template <typename T> struct PointAttribute : GenericAttribute<T> {
        PointAttribute(PointSet &pts, const T def = T());
        PointAttribute(const PointSet &pts, const T def = T());
        PointAttribute(PolyLine &m, const T def = T());
        PointAttribute(const PolyLine &m, const T def = T());
        PointAttribute(Surface &m, const T def = T());
        PointAttribute(const Surface &m, const T def = T());
        PointAttribute(Volume  &m, const T def = T());
        PointAttribute(const Volume  &m, const T def = T());

        PointAttribute(std::string name, PointSetAttributes &attributes, PointSet &ps,  const T def = T());
        PointAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def = T());
        PointAttribute(std::string name, SurfaceAttributes  &attributes, Surface  &m,   const T def = T());
        PointAttribute(std::string name, VolumeAttributes   &attributes, Volume   &m,   const T def = T());
        virtual AttributeBase::TYPE kind() const { return AttributeBase::POINTS; }
    };

    template <typename T> struct EdgeAttribute : GenericAttribute<T> {
        EdgeAttribute(PolyLine &seg, const T def = T());
        EdgeAttribute(const PolyLine &seg, const T def = T());
        EdgeAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def = T());
        virtual AttributeBase::TYPE kind() const { return AttributeBase::EDGES; }
    };

    template <typename T> struct FacetAttribute : GenericAttribute<T> {
        FacetAttribute(Surface &m, const T def = T());
        FacetAttribute(const Surface &m, const T def = T());
        FacetAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T());
        virtual AttributeBase::TYPE kind() const { return AttributeBase::FACETS; }
    };

    template <typename T> struct CornerAttribute : GenericAttribute<T> {
        CornerAttribute(Surface &m, const T def = T());
        CornerAttribute(const Surface &m, const T def = T());
        CornerAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T());
        virtual AttributeBase::TYPE kind() const { return AttributeBase::CORNERS; }
    };

    template <typename T> struct CellAttribute : GenericAttribute<T> {
        CellAttribute(Volume &m, const T def = T());
        CellAttribute(const Volume &m, const T def = T());
        CellAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T());
        virtual AttributeBase::TYPE kind() const { return AttributeBase::CELLS; }
    };

    template <typename T> struct CellFacetAttribute : GenericAttribute<T> {
        CellFacetAttribute(Volume &m, const T def = T());
        CellFacetAttribute(const Volume &m, const T def = T());
        CellFacetAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T());
        virtual AttributeBase::TYPE kind() const { return AttributeBase::CELLFACETS; }
    };

    template <typename T> struct CellCornerAttribute : GenericAttribute<T> {
        CellCornerAttribute(Volume &m, const T def = T());
        CellCornerAttribute(const Volume &m, const T def = T());
        CellCornerAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T());
        virtual AttributeBase::TYPE kind() const { return AttributeBase::CELLCORNERS; }
    };

    struct NamedContainer {
        std::string name;
        std::shared_ptr<ContainerBase> ptr;
    };

    struct NamedAttribute {
        std::string name;
        AttributeBase& attribute;
    };

    struct PointSetAttributes {
        PointSetAttributes() = default;
        PointSetAttributes(PointSetAttributes& p)        = default;
        PointSetAttributes(PointSetAttributes&& p)       = default;
        PointSetAttributes(const PointSetAttributes& p)  = default;
        PointSetAttributes& operator=(const PointSetAttributes& p)  = default;

        PointSetAttributes(std::vector<NamedContainer> list) : points(list) {}

        PointSetAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                um_assert(na.attribute.kind()==AttributeBase::POINTS);
                points.emplace_back(na.name, na.attribute.get_ptr());
            }
        }

        std::vector<NamedContainer> points;
    };

    struct PolyLineAttributes {
        PolyLineAttributes() = default;
        PolyLineAttributes(PolyLineAttributes& p)        = default;
        PolyLineAttributes(PolyLineAttributes&& p)       = default;
        PolyLineAttributes(const PolyLineAttributes& p)  = default;
        PolyLineAttributes& operator=(const PolyLineAttributes& p)  = default;

        PolyLineAttributes(std::vector<NamedContainer> points, std::vector<NamedContainer> edges) : points(points), edges(edges) {}

        PolyLineAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:  points.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::EDGES:    edges.emplace_back(na.name, na.attribute.get_ptr()); break;
                    default: um_assert(false);
                }
            }
        }

        std::vector<NamedContainer> points = {}, edges = {};
    };

    struct SurfaceAttributes {
        SurfaceAttributes() = default;
        SurfaceAttributes(SurfaceAttributes& p)        = default;
        SurfaceAttributes(SurfaceAttributes&& p)       = default;
        SurfaceAttributes(const SurfaceAttributes& p)  = default;
        SurfaceAttributes& operator=(const SurfaceAttributes& p)  = default;

        SurfaceAttributes(std::vector<NamedContainer> points, std::vector<NamedContainer> facets, std::vector<NamedContainer> corners) : points(points), facets(facets), corners(corners) {}

        SurfaceAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:   points.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::FACETS:   facets.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CORNERS: corners.emplace_back(na.name, na.attribute.get_ptr()); break;
                    default: um_assert(false);
                }
            }
        }

        std::vector<NamedContainer> points = {}, facets = {}, corners = {};
    };

    struct VolumeAttributes {
        VolumeAttributes() = default;
        VolumeAttributes(VolumeAttributes& p)        = default;
        VolumeAttributes(VolumeAttributes&& p)       = default;
        VolumeAttributes(const VolumeAttributes& p)  = default;
        VolumeAttributes& operator=(const VolumeAttributes& p)  = default;

        VolumeAttributes(std::vector<NamedContainer> points,
                         std::vector<NamedContainer> cells,
                         std::vector<NamedContainer> cell_facets,
                         std::vector<NamedContainer> cell_corners) : points(points),
                                                                     cells(cells),
                                                                     cell_facets(cell_facets),
                                                                     cell_corners(cell_corners) {
        }

        VolumeAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:   points.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CELLS:   cells.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CELLFACETS:   cell_facets.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CELLCORNERS:   cell_corners.emplace_back(na.name, na.attribute.get_ptr()); break;
                    default: um_assert(false);
                }
            }
        }

        std::vector<NamedContainer> points = {}, cells = {}, cell_facets = {}, cell_corners = {};
    };

}

#endif //__ATTRIBUTES_H__

