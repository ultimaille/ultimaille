namespace UM {

    inline auto PolyLine::iter_vertices() {
        struct iterator {
            Vertex v;
            void operator++() { ++v.id; }
            bool operator!=(const iterator& rhs) const { return v != rhs.v; }
            Vertex& operator*() { return v; }
        };
        struct wrapper {
            PolyLine& m;
            auto begin() { return iterator{ {m,0} }; }
            auto end() { return iterator{ {m,m.nverts()} }; }
        };
        return wrapper{ *this };
    }

    inline auto PolyLine::iter_edges() {
        struct iterator {
            Edge h;
            void operator++() {
                ++h.id;
                if (h.m.connected())
                    while (h.id < h.m.nedges() && !h.active()) ++h.id;
            }
            bool operator!=(const iterator& rhs) const { return h != rhs.h; }
            Edge& operator*() { return h; }
        };
        struct wrapper {
            PolyLine& m;
            auto begin() {
                iterator res{ {m,0} };
                if (m.connected())
                    if (m.nedges() > 0 && !res.h.active()) ++res;
                return res;
            }
            auto end() { return iterator{ {m,m.nedges()} }; }
        };
        return wrapper{ *this };
    }

    inline auto PolyLine::Vertex::iter_edges() {
        assert(m.connected());
        struct iterator {
            PolyLine::Edge he;
            void operator++() {
                const auto& c2c = he.m.conn->e2e;
                do {
                    he = c2c[he];
                    if (he == -1) return;
                } while (!he.active());
            }
            bool operator!=(const iterator& rhs) const { return he != rhs.he; }
            PolyLine::Edge& operator*() { return he; }
        };
        struct wrapper {
            Vertex v;
            iterator begin() { return { v.edge() }; }
            iterator end() { return { {v.m, -1} }; }
        };
        return wrapper{ *this };
    }

}

