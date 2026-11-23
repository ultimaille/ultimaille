namespace UM {

    inline auto Surface::iter_vertices() {
        struct iterator {
            Vertex v;
            iterator & operator++() { ++v.id; return *this; }
            bool operator!=(const iterator& rhs) const { return v != rhs.v; }
            Vertex& operator*() { return v; }
        };
        struct wrapper {
            Surface& m;
            iterator begin() { return {{m,0}}; }
            iterator end()   { return {{m,m.nverts()}}; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::iter_halfedges() {
        struct iterator {
            Halfedge h;
            const int ncorners; // need to freeze the end in the case the surface is changed in the meantime
            iterator & operator++() {
                ++h.id;
                if (h.m.connected())
                    while (h.id < ncorners && !h.active()) ++h.id;
                return *this;
            }
            bool operator!=(const iterator& rhs) const { return h != rhs.h; }
            Halfedge& operator*() { return h; }
        };
        struct wrapper {
            Surface& m;
            auto begin() { return ++iterator{{m, -1},           m.ncorners()}; } // ++(-1) is to select the first active
            auto end()   { return   iterator{{m, m.ncorners()}, m.ncorners()}; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::iter_facets() {
        struct iterator {
            Facet f;
            const int nfacets; // need to freeze the end in the case the surface is changed in the meantime
            iterator & operator++() {
                ++f.id;
                if (f.m.connected())
                    while (f.id < nfacets && !f.active()) ++f.id;
                return *this;
            }
            bool operator!=(const iterator& rhs) const { return f != rhs.f; }
            Facet& operator*() { return f; }
        };
        struct wrapper {
            Surface& m;
            auto begin() { return ++iterator{{m, -1},          m.nfacets()}; } // ++(-1) is to select the first active
            auto end()   { return   iterator{{m, m.nfacets()}, m.nfacets()}; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::Vertex::iter_halfedges() const {
        assert(m.connected());
        struct iterator {
            Surface::Halfedge he;
            void operator++() {
                const auto &c2c = he.m.conn->c2c;
                do {
                    he = c2c[he];
                    if (he == -1) return;
                } while (!he.active());
            }
            bool operator!=(const iterator& rhs) const { return he != rhs.he; }
            Surface::Halfedge& operator*() { return he; }
        };
        struct wrapper {
            Vertex v;
            iterator begin() { return { v.halfedge() }; }
            iterator end()   { return { {v.m, -1} }; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::Halfedge::iter_sector_halfedges() const {
        struct iterator {
            Surface::Halfedge h, seed;
            void operator++() {
                h = h.prev().opposite();
                if (h == seed) h = -1;
            }
            bool operator!=(const iterator& rhs) const { return h != rhs.h; }
            Surface::Halfedge& operator*() { return h; }
        };
        struct wrapper {
            Surface::Halfedge h;
            iterator begin() {
                const auto rewind_CCW = [](Surface::Halfedge h) -> Surface::Halfedge {
                    assert(h.active());
                    Halfedge cir = h;
                    do {
                        Halfedge opp = cir.opposite();
                        if (opp<=0) return cir;
                        cir = opp.next();
                    } while (cir != h);
                    return h;
                };
                Surface::Halfedge rwd = rewind_CCW(h);
                return {rwd, rwd};
            }
            iterator end() { return {{h.m, -1}, {h.m, -1}}; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::Facet::iter_halfedges() const {
        struct iterator {
            Surface::Halfedge he;
            void operator++() {
                auto f = he.facet();
                if (++he == f.halfedge()+f.size()) he = -1;
            }
            bool operator!=(const iterator& rhs) const { return he != rhs.he; }
            Surface::Halfedge& operator*() { return he; }
        };
        struct wrapper {
            Surface::Facet f;
            iterator begin() { return { f.halfedge() }; }
            iterator end()   { return { {f.m,-1} }; }
        };
        return wrapper{ *this };
    }

}

