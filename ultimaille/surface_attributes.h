#include <initializer_list>
#include "syntactic-sugar/assert.h"

namespace UM {
    struct SurfaceAttributes {
        SurfaceAttributes() = default;
        SurfaceAttributes(SurfaceAttributes& p)       = default;
        SurfaceAttributes(SurfaceAttributes&& p)      = default;
        SurfaceAttributes(const SurfaceAttributes& p) = default;
        SurfaceAttributes& operator=(const SurfaceAttributes& p) = default;

        SurfaceAttributes(AttributeMap points, AttributeMap facets, AttributeMap corners) : points(std::move(points)), facets(std::move(facets)), corners(std::move(corners)) {}

        SurfaceAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:   points[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::FACETS:   facets[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::CORNERS: corners[na.name] = na.attribute.get_ptr(); break;
                    default:
                                                 um_assert(false);
                }
            }
        }

        AttributeMap points  = {};
        AttributeMap facets  = {};
        AttributeMap corners = {};
    };

    template <typename T> struct Surface::FacetAttribute : GenericAttribute<T> {
        FacetAttribute(T def = T());
        FacetAttribute(Surface &m, T def = T());
        FacetAttribute(const Surface &m, T def = T());
        FacetAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def = T());
        bool bind(std::string name, SurfaceAttributes &attributes, Surface &m);
        virtual AttributeBase::TYPE kind() const { return AttributeBase::FACETS; }
    };

    template <typename T> struct Surface::CornerAttribute : GenericAttribute<T> {
        CornerAttribute(T def = T());
        CornerAttribute(Surface &m, T def = T());
        CornerAttribute(const Surface &m, T def = T());
        CornerAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def = T());
        bool bind(std::string name, SurfaceAttributes &attributes, Surface &m);
        virtual AttributeBase::TYPE kind() const { return AttributeBase::CORNERS; }
    };


    struct Surface::Connectivity {
        Surface& m;
        PointSet::Attribute<int>      v2c;    // vertex to corner map
        Surface::CornerAttribute<int> c2f;    // corner to facet map
        Surface::CornerAttribute<int> c2c;    // corner to corner (sharing the same vertex) map. Consecutive maps form an unordered linked list terminated by -1.
        Surface::FacetAttribute<bool> active; // facets to keep after compacting

        Connectivity(Surface& m);
        void init();
        Surface::Facet create_facet(std::initializer_list<int> verts);
        Surface::Facet create_facet(int* verts, int size);

        void change_from(Surface::Halfedge he, int new_vertex_id); // TODO move it to the toolbox
    };


    template <typename T> using FacetAttribute  = Surface::FacetAttribute<T>;
    template <typename T> using CornerAttribute = Surface::CornerAttribute<T>;

    template <typename T> PointSet::Attribute<T>::Attribute(const Surface &m, T def) : Attribute(m.points, def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(Surface  &m, T def) : Attribute(m.points, def) {}

    template <typename T> Surface::FacetAttribute<T>::FacetAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> Surface::FacetAttribute<T>::FacetAttribute(Surface &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
        m.attr_facets.push_back(this->ptr);
    }

    template <typename T> Surface::FacetAttribute<T>::FacetAttribute(const Surface &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
    }

    template <typename T> Surface::FacetAttribute<T>::FacetAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nfacets(), attributes.facets, m.attr_facets);
    }

    template <typename T> bool Surface::FacetAttribute<T>::bind(std::string name, SurfaceAttributes &attributes, Surface &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nfacets(), attributes.facets, m.attr_facets);
    }


    template <typename T> Surface::CornerAttribute<T>::CornerAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> Surface::CornerAttribute<T>::CornerAttribute(Surface &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
        m.attr_corners.push_back(this->ptr);
    }

    template <typename T> Surface::CornerAttribute<T>::CornerAttribute(const Surface &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
    }

    template <typename T> Surface::CornerAttribute<T>::CornerAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.ncorners(), attributes.corners, m.attr_corners);
    }

    template <typename T> bool Surface::CornerAttribute<T>::bind(std::string name, SurfaceAttributes &attributes, Surface &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.ncorners(), attributes.corners, m.attr_corners);
    }

    template <typename T> PointSet::Attribute<T>::Attribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nverts(), attributes.points, *m.points.attr);
    }

    template <typename T> bool PointSet::Attribute<T>::bind(std::string name, SurfaceAttributes &attributes, Surface &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr);
    }
}

