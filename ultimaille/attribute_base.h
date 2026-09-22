#ifndef __ATTRIBUTE_BASE_H__
#define __ATTRIBUTE_BASE_H__
#include <map>
#include <vector>
#include <memory>
#include "syntactic-sugar/assert.h"

namespace UM {
    struct ContainerBase {
        virtual ~ContainerBase() = default;
        virtual void resize(const int n) = 0;
        virtual void compress(const std::vector<int> &old2new) = 0;
        virtual void append_zeros(size_t n) = 0; // for merging attributes, at least in vtk.cpp
        virtual void prepend_zeros(size_t n) = 0;
        virtual void append(const ContainerBase& other) = 0;
    };

    struct AttributeBase {
        enum TYPE { GENERIC=-1, POINTS=0, EDGES=1, FACETS=2, CORNERS=3, CELLS=4, CELLFACETS=5, CELLCORNERS=6 };
        virtual ~AttributeBase() = default;
        virtual TYPE kind() const { return GENERIC; }
        virtual std::shared_ptr<ContainerBase> get_ptr() const = 0;
        bool bound() const { return get_ptr() != nullptr; }
    };

    template <typename T> struct AttributeContainer : ContainerBase {
        AttributeContainer(const int n, T def = T()) : data(n, def), default_value(def) {}
        void resize(const int n) { data.resize(n, default_value); }
        void compress(const std::vector<int> &old2new) { // NB: old2new is not a permutation!
            um_assert(old2new.size()==data.size());
            int cnt = 0;
            for (int i=0; i<static_cast<int>(old2new.size()); i++) {
                if (old2new[i]<0) continue;
                data[old2new[i]] = data[i];
                cnt++;
            }
            resize(cnt);
        }

        void append_zeros(std::size_t n) override { // TODO: not zeros, default value
            data.insert(data.end(), n, T{});
        }

        void prepend_zeros(std::size_t n) override {
            data.insert(data.begin(), n, T{});
        }

        void append(const ContainerBase& other) override {
            auto ptr = dynamic_cast<const AttributeContainer<T>*>(&other);
            um_assert(ptr != nullptr);
            data.insert(data.end(), ptr->data.begin(), ptr->data.end());
        }

        std::vector<T> data;
        T default_value;
    };

    template <typename T> struct GenericAttribute : AttributeBase {
        GenericAttribute(T default_value) : default_value(default_value), ptr(nullptr) {}
        GenericAttribute(T default_value, int size) : default_value(default_value), ptr(new AttributeContainer<T>(size, default_value)) {}

        GenericAttribute(const GenericAttribute<T>& rhs) = delete;
        GenericAttribute<T>& operator=(const GenericAttribute<T>& rhs) = delete;

              T& operator[](const int i)       { return ptr->data[i]; }
        const T& operator[](const int i) const { return ptr->data[i]; }
        void fill(T value) {
            if (ptr) std::fill(ptr->data.begin(), ptr->data.end(), value);
        }
        virtual std::shared_ptr<ContainerBase> get_ptr() const { return ptr; }

        T default_value; // necessary for the late binding
        std::shared_ptr<AttributeContainer<T> > ptr;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////

    template <> struct GenericAttribute<bool> : AttributeBase {
        GenericAttribute(bool default_value) : default_value(default_value), ptr(nullptr) {}
        GenericAttribute(bool default_value, int size) : default_value(default_value), ptr(new AttributeContainer<bool>(size, default_value)) {}

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

        bool default_value;
        std::shared_ptr<AttributeContainer<bool> > ptr;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////

    struct NamedAttribute {
        std::string name;
        AttributeBase& attribute;
    };

    using AttributeMap = std::map<std::string, std::shared_ptr<ContainerBase>>;

    inline void add_attribute(AttributeMap& map, const NamedAttribute& attribute, AttributeBase::TYPE expected) {
        um_assert(attribute.attribute.kind() == expected);
        map[attribute.name] = attribute.attribute.get_ptr();
    }
}

#endif //__ATTRIBUTE_BASE_H__

