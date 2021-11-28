#ifndef __ATTRIBUTES_H__
#define __ATTRIBUTES_H__
#include <vector>
#include <memory>
#include "polyline.h"
#include "surface.h"
#include "volume.h"
#include "pointset.h"

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

        std::shared_ptr<AttributeContainer<bool> > ptr;
    };

    template <typename T> void bind_attribute(GenericAttribute<T> *A, const std::string name, const int size, std::vector<NamedContainer> &containers, std::vector<std::weak_ptr<GenericAttributeContainer> > &callbacks, const T def = T()) {
        for (auto &pair : containers) {
            if (pair.first!=name) continue;
            A->ptr = std::dynamic_pointer_cast<AttributeContainer<T> >(pair.second);
            assert(A->ptr.get());
            A->ptr->default_value = def;
            //   callbacks.push_back(ptr); // TODO architectural choice: to bind or not to bind? At the moment the binding is done in mesh_io.cpp
            return;
        }
        A->ptr = std::make_shared<AttributeContainer<T> >(size, def);
        callbacks.push_back(A->ptr);
        containers.emplace_back(name, A->ptr);
    }

    template <typename T> struct PointAttribute : GenericAttribute<T> {
        PointAttribute(PointSet &pts, const T def = T()) : GenericAttribute<T>(pts.size(), def) {
            pts.attr.push_back(this->ptr);
        }

        PointAttribute(const PointSet &pts, const T def = T()) : GenericAttribute<T>(pts.size(), def) {
        }

        PointAttribute(Surface &m, const T def = T()) : PointAttribute(m.points, def) {}
        PointAttribute(const Surface &m, const T def = T()) : PointAttribute(m.points, def) {}
        PointAttribute(Volume  &m, const T def = T()) : PointAttribute(m.points, def) {}
        PointAttribute(const Volume  &m, const T def = T()) : PointAttribute(m.points, def) {}

        PointAttribute(std::string name, PointSetAttributes &attributes, PointSet &ps, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, ps.size(), attributes.points, ps.attr, def);
        }

        PointAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, seg.nverts(), attributes.points, seg.points.attr, def);
        }

        PointAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr, def);
        }

        PointAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr, def);
        }
    };

    template <typename T> struct SegmentAttribute : GenericAttribute<T> {
        SegmentAttribute(PolyLine &seg, const T def = T()) : GenericAttribute<T>(seg.nsegments(), def) {
            seg.attr.push_back(this->ptr);
        }

        SegmentAttribute(const PolyLine &seg, const T def = T()) : GenericAttribute<T>(seg.nsegments(), def) {
        }

        SegmentAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, seg.nsegments(), attributes.segments, seg.attr, def);
        }
    };

    template <typename T> struct FacetAttribute : GenericAttribute<T> {
        FacetAttribute(Surface &m, const T def = T()) : GenericAttribute<T>(m.nfacets(), def) {
            m.attr_facets.push_back(this->ptr);
        }

        FacetAttribute(const Surface &m, const T def = T()) : GenericAttribute<T>(m.nfacets(), def) {
        }

        FacetAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, m.nfacets(), attributes.facets, m.attr_facets, def);
        }
    };

    template <typename T> struct CornerAttribute : GenericAttribute<T> {
        CornerAttribute(Surface &m, const T def = T()) : GenericAttribute<T>(m.ncorners(), def) {
            m.attr_corners.push_back(this->ptr);
        }

        CornerAttribute(const Surface &m, const T def = T()) : GenericAttribute<T>(m.ncorners(), def) {
        }

        CornerAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, m.ncorners(), attributes.corners, m.attr_corners, def);
        }
    };

    template <typename T> struct CellAttribute : GenericAttribute<T> {
        CellAttribute(Volume &m, const T def = T()) : GenericAttribute<T>(m.ncells(), def) {
            m.attr_cells.push_back(this->ptr);
        }

        CellAttribute(const Volume &m, const T def = T()) : GenericAttribute<T>(m.ncells(), def) {
        }

        CellAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, m.ncells(), attributes.cells, m.attr_cells, def);
        }
    };

    template <typename T> struct CellFacetAttribute : GenericAttribute<T> {
        CellFacetAttribute(Volume &m, const T def = T()) : GenericAttribute<T>(m.nfacets(), def) {
            m.attr_facets.push_back(this->ptr);
        }

        CellFacetAttribute(const Volume &m, const T def = T()) : GenericAttribute<T>(m.nfacets(), def) {
        }

        CellFacetAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, m.nfacets(), attributes.cell_facets, m.attr_facets, def);
        }
    };

    template <typename T> struct CellCornerAttribute : GenericAttribute<T> {
        CellCornerAttribute(Volume &m, const T def = T()) : GenericAttribute<T>(m.ncorners(), def) {
            m.attr_corners.push_back(this->ptr);
        }

        CellCornerAttribute(const Volume &m, const T def = T()) : GenericAttribute<T>(m.ncorners(), def) {
        }

        CellCornerAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def = T()) : GenericAttribute<T>() {
            bind_attribute(this, name, m.ncorners(), attributes.cell_corners, m.attr_corners, def);
        }
    };
}

#endif //__ATTRIBUTES_H__

