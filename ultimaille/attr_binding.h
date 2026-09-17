#ifndef __ATTR_BINDING_H__
#define __ATTR_BINDING_H__
#include <vector>
#include <memory>

namespace UM {
    template <typename T> bool bind_attribute(GenericAttribute<T>* attribute, const std::string& name, int size, AttributeMap& containers, std::vector<std::weak_ptr<ContainerBase>>& callbacks) {
        auto it = containers.find(name);

        if (it != containers.end()) {
            attribute->ptr = std::dynamic_pointer_cast<AttributeContainer<T>>(it->second);
            um_assert(attribute->ptr != nullptr);
            attribute->ptr->default_value = attribute->default_value;
            return true;
        }

        attribute->ptr = std::make_shared<AttributeContainer<T>>(size, attribute->default_value);
        callbacks.push_back(attribute->ptr);
        containers.emplace(name, attribute->ptr);
        return false;
    }

    template <typename T> PointSet::Attribute<T>::Attribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(PointSet &pts, T def) : GenericAttribute<T>(def, pts.size()) {
        pts.attr->push_back(this->ptr);
    }
    template <typename T> PointSet::Attribute<T>::Attribute(const PointSet &pts, T def) : GenericAttribute<T>(def, pts.size()) {}

    template <typename T> PointSet::Attribute<T>::Attribute(PolyLine &m, T def) : Attribute(m.points, def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(Surface  &m, T def) : Attribute(m.points, def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(Volume   &m, T def) : Attribute(m.points, def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(const PolyLine &m, T def) : Attribute(m.points, def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(const Surface &m, T def) : Attribute(m.points, def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(const Volume  &m, T def) : Attribute(m.points, def) {}

    template <typename T> PointSet::Attribute<T>::Attribute(std::string name, PointSetAttributes &attributes, PointSet &ps, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, ps.size(), attributes.points, *ps.attr);
    }

    template <typename T> PointSet::Attribute<T>::Attribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, seg.nverts(), attributes.points, *seg.points.attr);
    }

    template <typename T> PointSet::Attribute<T>::Attribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nverts(), attributes.points, *m.points.attr);
    }

    template <typename T> PointSet::Attribute<T>::Attribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nverts(), attributes.points, *m.points.attr);
    }


    template <typename T> bool PointSet::Attribute<T>::bind(std::string name, PointSetAttributes &attributes, PointSet &ps) {
        um_assert(!this->bound());
        return bind_attribute(this, name, ps.size(), attributes.points, ps.attr);
    }

    template <typename T> bool PointSet::Attribute<T>::bind(std::string name, PolyLineAttributes &attributes, PolyLine &seg) {
        um_assert(!this->bound());
        return bind_attribute(this, name, seg.nverts(), attributes.points, seg.points.attr);
    }

    template <typename T> bool PointSet::Attribute<T>::bind(std::string name, SurfaceAttributes &attributes, Surface &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr);
    }

    template <typename T> bool PointSet::Attribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr);
    }


    template <typename T> PolyLine::Attribute<T>::Attribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> PolyLine::Attribute<T>::Attribute(PolyLine &seg, T def) : GenericAttribute<T>(def, seg.nedges()) {
        seg.attr.push_back(this->ptr);
    }

    template <typename T> PolyLine::Attribute<T>::Attribute(const PolyLine &seg, T def) : GenericAttribute<T>(def, seg.nedges()) {
    }

    template <typename T> PolyLine::Attribute<T>::Attribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, seg.nedges(), attributes.edges, seg.attr);
    }

    template <typename T> bool PolyLine::Attribute<T>::bind(std::string name, PolyLineAttributes &attributes, PolyLine &seg) {
        um_assert(!this->bound());
        return bind_attribute(this, name, seg.nedges(), attributes.edges, seg.attr);
    }

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


    template <typename T> Volume::CellAttribute<T>::CellAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> Volume::CellAttribute<T>::CellAttribute(Volume &m, T def) : GenericAttribute<T>(def, m.ncells()) {
        m.attr_cells.push_back(this->ptr);
    }

    template <typename T> Volume::CellAttribute<T>::CellAttribute(const Volume &m, T def) : GenericAttribute<T>(def, m.ncells()) {
    }

    template <typename T> Volume::CellAttribute<T>::CellAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.ncells(), attributes.cells, m.attr_cells);
    }

    template <typename T> bool Volume::CellAttribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.ncells(), attributes.cells, m.attr_cells);
    }

    template <typename T> Volume::FacetAttribute<T>::FacetAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> Volume::FacetAttribute<T>::FacetAttribute(Volume &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
        m.attr_facets.push_back(this->ptr);
    }

    template <typename T> Volume::FacetAttribute<T>::FacetAttribute(const Volume &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
    }

    template <typename T> Volume::FacetAttribute<T>::FacetAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nfacets(), attributes.cell_facets, m.attr_facets);
    }

    template <typename T> bool Volume::FacetAttribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nfacets(), attributes.cell_facets, m.attr_facets);
    }


    template <typename T> Volume::CornerAttribute<T>::CornerAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> Volume::CornerAttribute<T>::CornerAttribute(Volume &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
        m.attr_corners.push_back(this->ptr);
    }

    template <typename T> Volume::CornerAttribute<T>::CornerAttribute(const Volume &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
    }

    template <typename T> Volume::CornerAttribute<T>::CornerAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.ncorners(), attributes.cell_corners, m.attr_corners);
    }

    template <typename T> bool Volume::CornerAttribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.ncorners(), attributes.cell_corners, m.attr_corners);
    }
}

#endif //__ATTR_BINDING_H__

