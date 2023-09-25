#ifndef __ATTR_BINDING_H__
#define __ATTR_BINDING_H__
#include <vector>
#include <memory>

namespace UM {
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

    template <typename T> PointAttribute<T>::PointAttribute(PointSet &pts, const T def) : GenericAttribute<T>(pts.size(), def) {
        pts.attr.push_back(this->ptr);
    }
    template <typename T> PointAttribute<T>::PointAttribute(const PointSet &pts, const T def) : GenericAttribute<T>(pts.size(), def) {}

    template <typename T> PointAttribute<T>::PointAttribute(PolyLine &m, const T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(Surface  &m, const T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(Volume   &m, const T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(const PolyLine &m, const T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(const Surface &m, const T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(const Volume  &m, const T def) : PointAttribute(m.points, def) {}

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, PointSetAttributes &attributes, PointSet &ps, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, ps.size(), attributes.points, ps.attr, def);
    }

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, seg.nverts(), attributes.points, seg.points.attr, def);
    }

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr, def);
    }

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr, def);
    }




    template <typename T> EdgeAttribute<T>::EdgeAttribute(PolyLine &seg, const T def) : GenericAttribute<T>(seg.nedges(), def) {
        seg.attr.push_back(this->ptr);
    }

    template <typename T> EdgeAttribute<T>::EdgeAttribute(const PolyLine &seg, const T def) : GenericAttribute<T>(seg.nedges(), def) {
    }

    template <typename T> EdgeAttribute<T>::EdgeAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, seg.nedges(), attributes.edges, seg.attr, def);
    }




    template <typename T> FacetAttribute<T>::FacetAttribute(Surface &m, const T def) : GenericAttribute<T>(m.nfacets(), def) {
        m.attr_facets.push_back(this->ptr);
    }

    template <typename T> FacetAttribute<T>::FacetAttribute(const Surface &m, const T def) : GenericAttribute<T>(m.nfacets(), def) {
    }

    template <typename T> FacetAttribute<T>::FacetAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, m.nfacets(), attributes.facets, m.attr_facets, def);
    }




    template <typename T> CornerAttribute<T>::CornerAttribute(Surface &m, const T def) : GenericAttribute<T>(m.ncorners(), def) {
        m.attr_corners.push_back(this->ptr);
    }

    template <typename T> CornerAttribute<T>::CornerAttribute(const Surface &m, const T def) : GenericAttribute<T>(m.ncorners(), def) {
    }

    template <typename T> CornerAttribute<T>::CornerAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, m.ncorners(), attributes.corners, m.attr_corners, def);
    }




    template <typename T> CellAttribute<T>::CellAttribute(Volume &m, const T def) : GenericAttribute<T>(m.ncells(), def) {
        m.attr_cells.push_back(this->ptr);
    }

    template <typename T> CellAttribute<T>::CellAttribute(const Volume &m, const T def) : GenericAttribute<T>(m.ncells(), def) {
    }

    template <typename T> CellAttribute<T>::CellAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, m.ncells(), attributes.cells, m.attr_cells, def);
    }




    template <typename T> CellFacetAttribute<T>::CellFacetAttribute(Volume &m, const T def) : GenericAttribute<T>(m.nfacets(), def) {
        m.attr_facets.push_back(this->ptr);
    }

    template <typename T> CellFacetAttribute<T>::CellFacetAttribute(const Volume &m, const T def) : GenericAttribute<T>(m.nfacets(), def) {
    }

    template <typename T> CellFacetAttribute<T>::CellFacetAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, m.nfacets(), attributes.cell_facets, m.attr_facets, def);
    }




    template <typename T> CellCornerAttribute<T>::CellCornerAttribute(Volume &m, const T def) : GenericAttribute<T>(m.ncorners(), def) {
        m.attr_corners.push_back(this->ptr);
    }

    template <typename T> CellCornerAttribute<T>::CellCornerAttribute(const Volume &m, const T def) : GenericAttribute<T>(m.ncorners(), def) {
    }

    template <typename T> CellCornerAttribute<T>::CellCornerAttribute(std::string name, VolumeAttributes &attributes, Volume &m, const T def) : GenericAttribute<T>() {
        bind_attribute(this, name, m.ncorners(), attributes.cell_corners, m.attr_corners, def);
    }
}

#endif //__ATTR_BINDING_H__

