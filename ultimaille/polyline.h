#ifndef __POLYLINE_H__
#define __POLYLINE_H__

#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "pointset.h"
#include "attributes.h"
#include "syntactic-sugar/assert.h"
#include "primitive_geometry.h"

namespace UM {
    struct GenericAttributeContainer;

    struct PolyLine {
        PointSet points{};
        std::vector<int> edges{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr{};

        int nverts() const;
        int nedges() const;

        void compress_attrs(const std::vector<bool> &edges_to_kill);
        void delete_vertices(const std::vector<bool> &to_kill);
        void delete_edges(const std::vector<bool> &to_kill);
        int create_edges(const int n);
        void resize_attrs();

        int  vert(const int s, const int lv) const;
        int &vert(const int s, const int lv)      ;

        virtual void clear() {
            points   = {};
            attr     = {};
            edges = {};
            disconnect();
        }

        void delete_isolated_vertices();

        PolyLine() {}
        PolyLine(const PolyLine& m) {
            um_assert(!m.points.size() && !m.edges.size());
        }
        PolyLine& operator=(const PolyLine& m) {
            clear();
            um_assert(!m.points.size() && !m.edges.size());
            return *this;
        }

        struct Edge;
        struct Connectivity {
            PolyLine& m;
            PointAttribute<int>     v2e;    // vertex to edge
            EdgeAttribute<int>      e2e;    // edge to next edge sharing the same origin (unordered)
            EdgeAttribute<bool>     active; // edges to keep after compacting

            Connectivity(PolyLine& m);
            void init();
            Edge create_edge(int v0, int v1);
        };

        std::unique_ptr<Connectivity> conn = {};
        inline bool connected() const { return conn != nullptr; }

        void connect();
        void disconnect();
        void compact(bool delete_isolated_vertices = true);

        struct Primitive {
            Primitive(PolyLine& m, int id);
            Primitive(Primitive& p) = default;
            Primitive(Primitive&& p) = default;

            Primitive& operator=(Primitive& p);
            Primitive& operator=(int i);

            operator int() const;
            operator int& ();

        protected:
            friend struct PolyLine;
            PolyLine& m;
            int id;
        };


        struct Vertex : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Vertex(Vertex& v) = default;
            Vertex(Vertex&& v) = default;
            Vertex& operator=(Vertex& v);

            vec3  pos() const;
            vec3& pos();
            Edge edge();
            bool active();
            auto iter_edges();
        };

        struct Edge : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Edge(Edge& he) = default;
            Edge(Edge&& he) = default;
            Edge& operator=(Edge& he);

            bool active();
            Edge opposite();
            Vertex from() const;
            Vertex to() const;

            [[deprecated]] inline Segment3 geom();
            operator Segment3() const;

            friend struct Vertex;
        };

        Vertex vertex(int id) { return Vertex(*this, id); }
        Edge edge(int id) { return Edge(*this, id); }

        auto iter_vertices();
        auto iter_edges();
    };



/*   _ _                 _
    (_) |               | |
     _| |_ ___ _ __ __ _| |_ ___  _ __ ___
    | | __/ _ \ '__/ _` | __/ _ \| '__/ __|
    | | ||  __/ | | (_| | || (_) | |  \__ \
    |_|\__\___|_|  \__,_|\__\___/|_|  |___/  */



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




/*  ______     _           _ _   _
    | ___ \   (_)         (_) | (_)
    | |_/ / __ _ _ __ ___  _| |_ ___   _____  ___
    |  __/ '__| | '_ ` _ \| | __| \ \ / / _ \/ __|
    | |  | |  | | | | | | | | |_| |\ V /  __/\__ \
    \_|  |_|  |_|_| |_| |_|_|\__|_| \_/ \___||___/  */

inline PolyLine::Primitive::Primitive(PolyLine& m, int id) : m(m), id(id) {}

inline PolyLine::Primitive& PolyLine::Primitive::operator=(PolyLine::Primitive& p) {
    assert(&m == &p.m);
    id = p.id;
    return *this;
}

inline PolyLine::Primitive& PolyLine::Primitive::operator=(int i) {
    id = i;
    return *this;
}

inline PolyLine::Primitive::operator int() const { return id; }
inline PolyLine::Primitive::operator int& () { return id; }

////////////////////////////////////////////////////////////////////////////////////////////////////////

inline PolyLine::Vertex& PolyLine::Vertex::operator=(PolyLine::Vertex& v) {
    Primitive::operator=(v);
    return *this;
}

inline vec3  PolyLine::Vertex::pos() const { return m.points[id]; }
inline vec3& PolyLine::Vertex::pos() { return m.points[id]; }

inline PolyLine::Edge PolyLine::Vertex::edge() {
    assert(m.connected());
    PolyLine::Edge res{m, m.conn->v2e[id]};
    while (res >= 0 && !res.active()) {
        res = m.conn->e2e[res];
    }
    return res;
}

inline bool PolyLine::Vertex::active() {
    return id >= 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

inline PolyLine::Edge& PolyLine::Edge::operator=(PolyLine::Edge& he) {
    Primitive::operator=(he);
    return *this;
}

inline bool PolyLine::Edge::active() {
    return id >= 0 && m.conn->active[id];
}


inline PolyLine::Edge PolyLine::Edge::opposite() {
    assert(m.connected());
    assert(active());
    Edge result = { m, -1 }; // not found
    for (auto candidate : to().iter_edges()) {
        um_assert(candidate.from() == to());
        if (candidate.to() == from()) {
            if (result == -1) result = candidate;
            else return { m, -1 }; // found more than one
        }
    }
    return result;
}

inline PolyLine::Vertex PolyLine::Edge::from() const {
    return { m, m.vert(id, 0) };
}

inline PolyLine::Vertex PolyLine::Edge::to() const {
    return { m, m.vert(id, 1) };
}

inline Segment3 PolyLine::Edge::geom() {
    return Segment3(from().pos(), to().pos());
}

inline PolyLine::Edge::operator Segment3() const {
    return {from().pos(), to().pos()};
}



}

#endif //__POLYLINE_H__
