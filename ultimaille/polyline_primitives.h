#include "primitive_geometry.h"

namespace UM {
    struct PolyLine::Connectivity {
        PolyLine& m;
        PointSet::Attribute<int>  v2e;    // vertex to edge
        PolyLine::Attribute<int>  e2e;    // edge to next edge sharing the same origin (unordered)
        PolyLine::Attribute<bool> active; // edges to keep after compacting

        Connectivity(PolyLine& m);
        void init();
        Edge create_edge(int v0, int v1);
    };

    inline bool PolyLine::connected() const { return conn != nullptr; }

    struct PolyLine::Primitive {
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

    struct PolyLine::Vertex : Primitive {
        using Primitive::Primitive;
        using Primitive::operator=;
        Vertex(Vertex& v) = default;
        Vertex(Vertex&& v) = default;
        Vertex& operator=(Vertex& v);

        inline operator vec3&();
        inline operator vec3&() const;
        vec3  pos() const;
        vec3& pos();
        Edge edge();
        bool active();
        auto iter_edges();
    };

    struct PolyLine::Edge : Primitive {
        using Primitive::Primitive;
        using Primitive::operator=;
        Edge(Edge& he) = default;
        Edge(Edge&& he) = default;
        Edge& operator=(Edge& he);

        bool active();
        Edge opposite();
        Vertex from() const;
        Vertex to() const;
        Vertex vertex(int lv) const;

        operator Segment3() const;

        friend struct Vertex;
    };

    inline PolyLine::Vertex PolyLine::vertex(int id) { return Vertex(*this, id); }
    inline PolyLine::Edge   PolyLine::edge  (int id) { return   Edge(*this, id); }

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
    inline PolyLine::Primitive::operator int&()      { return id; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline PolyLine::Vertex& PolyLine::Vertex::operator=(PolyLine::Vertex& v) {
        Primitive::operator=(v);
        return *this;
    }
    inline PolyLine::Vertex::operator vec3&() {
        return { m.points[id] };
    }

    inline PolyLine::Vertex::operator vec3&() const {
        return { m.points[id] };
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

    inline PolyLine::Vertex PolyLine::Edge::from() const {
        return { m, m.vert(id, 0) };
    }

    inline PolyLine::Vertex PolyLine::Edge::to() const {
        return { m, m.vert(id, 1) };
    }

    inline PolyLine::Vertex PolyLine::Edge::vertex(int lv) const {
        assert(lv == 0 || lv == 1);
        return { m, m.vert(id, lv) };
    }

    inline PolyLine::Edge::operator Segment3() const {
        return {from(), to()};
    }
}

