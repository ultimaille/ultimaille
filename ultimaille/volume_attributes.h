template <typename T> struct Volume::CellAttribute : GenericAttribute<T> {
    CellAttribute(T def = T());
    CellAttribute(Volume &m, T def = T());
    CellAttribute(const Volume &m, T def = T());
    CellAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def = T());
    bool bind(std::string name, VolumeAttributes &attributes, Volume &m);
    virtual AttributeBase::TYPE kind() const { return AttributeBase::CELLS; }
};

template <typename T> struct Volume::FacetAttribute : GenericAttribute<T> {
    FacetAttribute(T def = T());
    FacetAttribute(Volume &m, T def = T());
    FacetAttribute(const Volume &m, T def = T());
    FacetAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def = T());
    bool bind(std::string name, VolumeAttributes &attributes, Volume &m);
    virtual AttributeBase::TYPE kind() const { return AttributeBase::CELLFACETS; }
};

template <typename T> struct Volume::CornerAttribute : GenericAttribute<T> {
    CornerAttribute(T def = T());
    CornerAttribute(Volume &m, T def = T());
    CornerAttribute(const Volume &m, T def = T());
    CornerAttribute(std::string name, VolumeAttributes &attributes, Volume &m, T def = T());
    bool bind(std::string name, VolumeAttributes &attributes, Volume &m);
    virtual AttributeBase::TYPE kind() const { return AttributeBase::CELLCORNERS; }
};

template <typename T> using CellAttribute       = Volume::CellAttribute<T>;
template <typename T> using CellFacetAttribute  = Volume::FacetAttribute<T>;
template <typename T> using CellCornerAttribute = Volume::CornerAttribute<T>;

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct VolumeAttributes {
    VolumeAttributes() = default;
    VolumeAttributes(VolumeAttributes& p)       = default;
    VolumeAttributes(VolumeAttributes&& p)      = default;
    VolumeAttributes(const VolumeAttributes& p) = default;
    VolumeAttributes& operator=(const VolumeAttributes& p) = default;

    VolumeAttributes(
            AttributeMap points,
            AttributeMap cells,
            AttributeMap cell_facets,
            AttributeMap cell_corners) :
        points(std::move(points)),
        cells(std::move(cells)),
        cell_facets(std::move(cell_facets)),
        cell_corners(std::move(cell_corners)) {}

    VolumeAttributes(std::initializer_list<NamedAttribute> list) {
        for (const auto& na : list) {
            switch (na.attribute.kind()) {
                case AttributeBase::POINTS:            points[na.name] = na.attribute.get_ptr(); break;
                case AttributeBase::CELLS:              cells[na.name] = na.attribute.get_ptr(); break;
                case AttributeBase::CELLFACETS:   cell_facets[na.name] = na.attribute.get_ptr(); break;
                case AttributeBase::CELLCORNERS: cell_corners[na.name] = na.attribute.get_ptr(); break;
                default:
                                                 um_assert(false);
            }
        }
    }

    AttributeMap points       = {};
    AttributeMap cells        = {};
    AttributeMap cell_facets  = {};
    AttributeMap cell_corners = {};
};

template <typename T> PointSet::Attribute<T>::Attribute(const Volume  &m, T def) : Attribute(m.points, def) {}
template <typename T> PointSet::Attribute<T>::Attribute(Volume   &m, T def) : Attribute(m.points, def) {}

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

template <typename T> PointSet::Attribute<T>::Attribute(std::string name, VolumeAttributes &attributes, Volume &m, T def) : GenericAttribute<T>(def) {
    bind_attribute(this, name, m.nverts(), attributes.points, *m.points.attr);
}

template <typename T> bool PointSet::Attribute<T>::bind(std::string name, VolumeAttributes &attributes, Volume &m) {
    um_assert(!this->bound());
    return bind_attribute(this, name, m.nverts(), attributes.points, *m.points.attr);
}

