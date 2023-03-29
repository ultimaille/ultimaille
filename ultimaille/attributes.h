#ifndef __ATTRIBUTES_H__
#define __ATTRIBUTES_H__
#include <vector>
#include <memory>
#include <cassert>
#include "polyline.h"
#include "pointset.h"
#include "surface.h"
#include "volume.h"

namespace UM {
    struct GenericAttributeContainer {
        virtual void resize(const int n) = 0;
        virtual void compress(const std::vector<int> &old2new) = 0;
        virtual ~GenericAttributeContainer() = default;
    };

    template <typename T> struct AttributeContainer : GenericAttributeContainer {
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

    typedef std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > NamedContainer;
    struct PointSetAttributes {
        std::vector<NamedContainer> points;
    };
    struct PolyLineAttributes {
        std::vector<NamedContainer> points, segments;
    };
    struct SurfaceAttributes {
        std::vector<NamedContainer> points, facets, corners;
    };
    struct VolumeAttributes {
        std::vector<NamedContainer> points, cells, cell_facets, cell_corners;
    };

    template <typename T> struct GenericAttribute {
        GenericAttribute() : ptr(nullptr) {}
        GenericAttribute(int size, const T def = T()) : ptr(new AttributeContainer<T>(size, def)) {}
        GenericAttribute(std::shared_ptr<AttributeContainer<T> > p) : ptr(p) {}
        GenericAttribute(const GenericAttribute<T>& rhs) = delete;
        GenericAttribute<bool>& operator=(const GenericAttribute<T>& rhs) = delete;
              T& operator[](const int i)       { return ptr->data[i]; }
        const T& operator[](const int i) const { return ptr->data[i]; }
        void fill(T value) {
            if (ptr) std::fill(ptr->data.begin(), ptr->data.end(), value);
        }
        std::shared_ptr<AttributeContainer<T> > ptr;
    };

    template <> struct GenericAttribute<bool> {
        GenericAttribute() : ptr(nullptr) {}
        GenericAttribute(int size, const bool def = false) : ptr(new AttributeContainer<bool>(size, def)) {}
        GenericAttribute(std::shared_ptr<AttributeContainer<bool> > p) : ptr(p) {}
        GenericAttribute(const GenericAttribute<bool>& rhs) = delete;
        GenericAttribute<bool>& operator=(const GenericAttribute<bool>& rhs) = delete;

        struct ConstBoolAttributeAccessor {
            ConstBoolAttributeAccessor(const GenericAttribute<bool>& attribute, int index) : attribute_(&attribute), index_(index) {}
            operator bool() const {
                return attribute_->ptr->data[index_];
            }
            const GenericAttribute<bool>* attribute_;
            const int index_;
        };

        struct BoolAttributeAccessor {
            BoolAttributeAccessor(GenericAttribute<bool>& attribute, int index) : attribute_(&attribute), index_(index) { }
            operator bool() const {
                return attribute_->ptr->data[index_];
            }

            BoolAttributeAccessor(const BoolAttributeAccessor& rhs) {
                attribute_ = rhs.attribute_;
                index_ = rhs.index_;
            }

            BoolAttributeAccessor& operator=(bool x) {
                attribute_->ptr->data[index_] = x;
                return *this;
            }

            BoolAttributeAccessor& operator=(const BoolAttributeAccessor& rhs) {
                if (&rhs != this)
                    attribute_->ptr->data[index_] = rhs.attribute_->ptr->data[rhs.index_];
                return *this;
            }

            BoolAttributeAccessor& operator=(const ConstBoolAttributeAccessor& rhs) {
                attribute_->ptr->data[index_] = rhs.attribute_->ptr->data[rhs.index_];
                return *this;
            }

            GenericAttribute<bool>* attribute_;
            int index_;
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

        PointAttribute(std::string name, PointSetAttributes &attributes, PointSet &ps, const T def = T());
        PointAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def = T());
        PointAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T());
        PointAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T());
    };

    template <typename T> struct SegmentAttribute : GenericAttribute<T> {
        SegmentAttribute(PolyLine &seg, const T def = T());
        SegmentAttribute(const PolyLine &seg, const T def = T());
        SegmentAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def = T());
    };

    template <typename T> struct FacetAttribute : GenericAttribute<T> {
        FacetAttribute(Surface &m, const T def = T());
        FacetAttribute(const Surface &m, const T def = T());
        FacetAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T());
    };

    template <typename T> struct CornerAttribute : GenericAttribute<T> {
        CornerAttribute(Surface &m, const T def = T());
        CornerAttribute(const Surface &m, const T def = T());
        CornerAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T());
    };

    template <typename T> struct CellAttribute : GenericAttribute<T> {
        CellAttribute(Volume &m, const T def = T());
        CellAttribute(const Volume &m, const T def = T());
        CellAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T());
    };

    template <typename T> struct CellFacetAttribute : GenericAttribute<T> {
        CellFacetAttribute(Volume &m, const T def = T());
        CellFacetAttribute(const Volume &m, const T def = T());
        CellFacetAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T());
    };

    template <typename T> struct CellCornerAttribute : GenericAttribute<T> {
        CellCornerAttribute(Volume &m, const T def = T());
        CellCornerAttribute(const Volume &m, const T def = T());
        CellCornerAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T());
    };
}

#endif //__ATTRIBUTES_H__

