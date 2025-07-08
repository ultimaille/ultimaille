#ifndef __ATTR_BINDING_H__
#define __ATTR_BINDING_H__
#include <vector>
#include <memory>

namespace UM {
    template <typename T> bool bind_attribute(GenericAttribute<T> *A, const std::string name, const int size, std::vector<NamedContainer> &containers, std::vector<std::weak_ptr<ContainerBase> > &callbacks) {
        for (auto &pair : containers) {
            if (pair.name!=name) continue;
            A->ptr = std::dynamic_pointer_cast<AttributeContainer<T> >(pair.ptr);
            assert(A->ptr.get());
            A->ptr->default_value = A->default_value;
            //   callbacks.push_back(ptr); // TODO architectural choice: to bind or not to bind? At the moment the binding is done in mesh_io.cpp
            return true;
        }
        A->ptr = std::make_shared<AttributeContainer<T> >(size, A->default_value);
        callbacks.push_back(A->ptr);
        containers.emplace_back(name, A->ptr);
        return false;
    }

    template <typename T> PointAttribute<T>::PointAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> PointAttribute<T>::PointAttribute(PointSet &pts, T def) : GenericAttribute<T>(def, pts.size()) {
        pts.attr.push_back(this->ptr);
    }
    template <typename T> PointAttribute<T>::PointAttribute(const PointSet &pts, T def) : GenericAttribute<T>(def, pts.size()) {}

    template <typename T> PointAttribute<T>::PointAttribute(PolyLine &m, T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(Surface  &m, T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(Volume   &m, T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(const PolyLine &m, T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(const Surface &m, T def) : PointAttribute(m.points, def) {}
    template <typename T> PointAttribute<T>::PointAttribute(const Volume  &m, T def) : PointAttribute(m.points, def) {}

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, PointSetAttributes &attributes, PointSet &ps, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, ps.size(), attributes.points, ps.attr);
    }

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, seg.nverts(), attributes.points, seg.points.attr);
    }

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr);
    }

    template <typename T> PointAttribute<T>::PointAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr);
    }


    template <typename T> bool PointAttribute<T>::bind(std::string name, PointSetAttributes &attributes, PointSet &ps) {
        um_assert(!this->bound());
        return bind_attribute(this, name, ps.size(), attributes.points, ps.attr);
    }

    template <typename T> bool PointAttribute<T>::bind(std::string name, PolyLineAttributes &attributes, PolyLine &seg) {
        um_assert(!this->bound());
        return bind_attribute(this, name, seg.nverts(), attributes.points, seg.points.attr);
    }

    template <typename T> bool PointAttribute<T>::bind(std::string name, SurfaceAttributes &attributes, Surface &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr);
    }

    template <typename T> bool PointAttribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nverts(), attributes.points, m.points.attr);
    }



    template <typename T> EdgeAttribute<T>::EdgeAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> EdgeAttribute<T>::EdgeAttribute(PolyLine &seg, T def) : GenericAttribute<T>(def, seg.nedges()) {
        seg.attr.push_back(this->ptr);
    }

    template <typename T> EdgeAttribute<T>::EdgeAttribute(const PolyLine &seg, T def) : GenericAttribute<T>(def, seg.nedges()) {
    }

    template <typename T> EdgeAttribute<T>::EdgeAttribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, seg.nedges(), attributes.edges, seg.attr);
    }

    template <typename T> bool EdgeAttribute<T>::bind(std::string name, PolyLineAttributes &attributes, PolyLine &seg) {
        um_assert(!this->bound());
        return bind_attribute(this, name, seg.nedges(), attributes.edges, seg.attr);
    }





    template <typename T> FacetAttribute<T>::FacetAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> FacetAttribute<T>::FacetAttribute(Surface &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
        m.attr_facets.push_back(this->ptr);
    }

    template <typename T> FacetAttribute<T>::FacetAttribute(const Surface &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
    }

    template <typename T> FacetAttribute<T>::FacetAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nfacets(), attributes.facets, m.attr_facets);
    }

    template <typename T> bool FacetAttribute<T>::bind(std::string name, SurfaceAttributes &attributes, Surface &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nfacets(), attributes.facets, m.attr_facets);
    }




    template <typename T> CornerAttribute<T>::CornerAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> CornerAttribute<T>::CornerAttribute(Surface &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
        m.attr_corners.push_back(this->ptr);
    }

    template <typename T> CornerAttribute<T>::CornerAttribute(const Surface &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
    }

    template <typename T> CornerAttribute<T>::CornerAttribute(std::string name, SurfaceAttributes &attributes, Surface &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.ncorners(), attributes.corners, m.attr_corners);
    }

    template <typename T> bool CornerAttribute<T>::bind(std::string name, SurfaceAttributes &attributes, Surface &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.ncorners(), attributes.corners, m.attr_corners);
    }




    template <typename T> CellAttribute<T>::CellAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> CellAttribute<T>::CellAttribute(Volume &m, T def) : GenericAttribute<T>(def, m.ncells()) {
        m.attr_cells.push_back(this->ptr);
    }

    template <typename T> CellAttribute<T>::CellAttribute(const Volume &m, T def) : GenericAttribute<T>(def, m.ncells()) {
    }

    template <typename T> CellAttribute<T>::CellAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.ncells(), attributes.cells, m.attr_cells);
    }

    template <typename T> bool CellAttribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.ncells(), attributes.cells, m.attr_cells);
    }




    template <typename T> CellFacetAttribute<T>::CellFacetAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> CellFacetAttribute<T>::CellFacetAttribute(Volume &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
        m.attr_facets.push_back(this->ptr);
    }

    template <typename T> CellFacetAttribute<T>::CellFacetAttribute(const Volume &m, T def) : GenericAttribute<T>(def, m.nfacets()) {
    }

    template <typename T> CellFacetAttribute<T>::CellFacetAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nfacets(), attributes.cell_facets, m.attr_facets);
    }

    template <typename T> bool CellFacetAttribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nfacets(), attributes.cell_facets, m.attr_facets);
    }




    template <typename T> CellCornerAttribute<T>::CellCornerAttribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> CellCornerAttribute<T>::CellCornerAttribute(Volume &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
        m.attr_corners.push_back(this->ptr);
    }

    template <typename T> CellCornerAttribute<T>::CellCornerAttribute(const Volume &m, T def) : GenericAttribute<T>(def, m.ncorners()) {
    }

    template <typename T> CellCornerAttribute<T>::CellCornerAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.ncorners(), attributes.cell_corners, m.attr_corners);
    }

    template <typename T> bool CellCornerAttribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.ncorners(), attributes.cell_corners, m.attr_corners);
    }
}

#endif //__ATTR_BINDING_H__

