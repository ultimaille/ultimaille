template <typename P, typename M>
auto custom_iterator(M &m, int from, int to) {
    struct iterator {
        P data;
        P& operator*()    { return data; }
        void operator++() { ++data; }
        bool operator!=(iterator& rhs) { return data != rhs.data; }
    };
    struct wrapper {
        M& m;
        int from, to;
        iterator begin() { return {{ m, from }}; }
        iterator end()   { return {{ m, to   }}; }
    };
    return wrapper{ m, from, to };
}

inline auto Volume::iter_vertices() const {
    return custom_iterator<Vertex>(*this, 0, nverts());
}

inline auto Volume::iter_corners() const {
    return custom_iterator<Corner>(*this, 0, ncorners());
}

inline auto Volume::iter_halfedges() const {
    return custom_iterator<Halfedge>(*this, 0, nhalfedges());
}

inline auto Volume::iter_facets() const {
    return custom_iterator<Facet>(*this, 0, nfacets());
}

inline auto Volume::iter_cells() const {
    return custom_iterator<Cell>(*this, 0, ncells());
}

inline auto Volume::Facet::iter_halfedges() {
    return custom_iterator<Halfedge>(m, halfedge(0), halfedge(0) + size());
}

inline auto Volume::Cell::iter_halfedges() {
    return custom_iterator<Halfedge>(m, halfedge(0), halfedge(0) + nhalfedges());
}

inline auto Volume::Cell::iter_facets() {
    return custom_iterator<Facet>(m, facet(0), facet(0) + nfacets());
}

inline auto Volume::Cell::iter_corners() {
    return custom_iterator<Corner>(m, corner(0), corner(0) + ncorners());
}

