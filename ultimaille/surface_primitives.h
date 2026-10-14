#include "primitive_geometry.h"

namespace UM {

    struct Surface::Primitive {
        Primitive(Surface& m, int id);
        Primitive(Primitive& p)        = default;
        Primitive(Primitive&& p)       = default;
        Primitive(const Primitive& p)  = default;

        Primitive& operator=(const Primitive&& p) noexcept;
        Primitive& operator=(const Primitive& p);
        Primitive& operator=(int i);

        operator int() const;
        operator int&();

        protected:
        friend struct Surface;
        Surface& m;
        int id;
    };

    struct Surface::Vertex : Primitive {
        using Primitive::Primitive;
        using Primitive::operator=;

        inline operator vec3&();
        inline operator vec3&() const;
        vec3  pos() const;
        vec3& pos();
        Halfedge halfedge() const;
        bool on_boundary() const;
        bool active() const;
        int id_in_facet(Facet f) const;

        auto iter_halfedges() const;
    };

    struct Surface::Halfedge : Primitive {
        using Primitive::Primitive;
        using Primitive::operator=;

        bool active() const;
        bool on_boundary() const;
        Facet facet() const;

        Halfedge next() const;
        Halfedge prev() const;
        Halfedge opposite() const;

        int id_in_facet() const;

        friend struct Vertex;
        Vertex from() const;
        Vertex to() const;

        inline operator vec3() const;
        inline operator Segment3() const;

        auto iter_sector_halfedges() const;
    };

    struct Surface::Facet : Primitive {
        using Primitive::Primitive;
        using Primitive::operator=;

        void deactivate();

        Vertex vertex(int lv) const;
        Halfedge halfedge(int lh = 0) const;
        int size() const;
        bool active() const;

        operator Triangle3() const;
        operator Quad3() const;
        operator Poly3() const;

        auto iter_halfedges() const;
    };


    inline Surface::Vertex     Surface::vertex(int id) { return   Vertex(*this, id); }
    inline Surface::Halfedge Surface::halfedge(int id) { return Halfedge(*this, id); }
    inline Surface::Facet       Surface::facet(int id) { return    Facet(*this, id); }


    ////////////////////////////////////////////////////////////////////////
    //      _____                            _   _       _ _         _    //
    //     / ____|                          | | (_)     (_) |       | |   //
    //    | |     ___  _ __  _ __   ___  ___| |_ ___   ___| |_ _   _| |   //
    //    | |    / _ \| '_ \| '_ \ / _ \/ __| __| \ \ / / | __| | | | |   //
    //    | |___| (_) | | | | | | |  __/ (__| |_| |\ V /| | |_| |_| |_|   //
    //     \_____\___/|_| |_|_| |_|\___|\___|\__|_| \_/ |_|\__|\__, (_)   //
    //                                                          __/ |     //
    //                                                         |___/      //
    ////////////////////////////////////////////////////////////////////////

    inline Surface::Primitive::Primitive(Surface& m, int id) : m(m), id(id) {}

    inline Surface::Primitive& Surface::Primitive::operator=(const Surface::Primitive&& p) noexcept {
        return Primitive::operator=(p);
    }

    inline Surface::Primitive& Surface::Primitive::operator=(const Surface::Primitive& p) {
        assert(&m == &p.m);
        id = p.id;
        return *this;
    }

    inline Surface::Primitive& Surface::Primitive::operator=(int i) {
        id = i;
        return *this;
    }

    inline Surface::Primitive::operator int() const { return id; }
    inline Surface::Primitive::operator int&()      { return id; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Surface::Vertex::operator vec3&() {
        return { m.points[id] };
    }

    inline Surface::Vertex::operator vec3&() const {
        return { m.points[id] };
    }

    inline vec3  Surface::Vertex::pos() const { return m.points[id]; }
    inline vec3& Surface::Vertex::pos()       { return m.points[id]; }

    inline Surface::Halfedge Surface::Vertex::halfedge() const {
        assert(m.connected());
        Surface::Halfedge res{m, m.conn->v2c[id]};
        while (res >= 0 && !res.active())
            res = m.conn->c2c[res];
        return res;
    }

    inline bool Surface::Vertex::active() const {
        return id >= 0;
    }

    inline int Surface::Vertex::id_in_facet(Facet f) const {
        for (int lv = 0; lv < f.size(); lv++)
            if (m.vert(f, lv) == id) return lv;
        assert(false);
        return -1;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline bool Surface::Halfedge::active() const {
        return id >= 0 && facet().active();
    }

    inline bool Surface::Halfedge::on_boundary() const {
        assert(m.connected());
        assert(active());

        int result = -1; // not found
        Halfedge cir = from().halfedge();
        while (cir != -1) {
            Halfedge candidate = cir.prev();
            if (cir.active()) {
                if (candidate.from() == to() && candidate.to() == from()) {
                    if (result == -1) result = candidate;
                    // found more than one (non manifold)
                    else
                        assert(false);
                }
                // the edge is non manifold (non-orientable)
                if (cir != *this && to() == cir.to())
                    assert(false);
            }
            do {
                cir = m.conn->c2c[cir];
            } while (cir != -1 && !cir.active());
        }

        // Found an opposite
        return result < 0;
    }

    inline Surface::Facet Surface::Halfedge::facet() const {
        assert(m.connected());
        return { m, m.conn->c2f[id] };
    }

    inline Surface::Halfedge Surface::Halfedge::next() const {
        auto f = facet();
        int lh = id - f.halfedge();
        int n = f.size();
        return { m, f.halfedge((lh + 1) % n) };
    }

    inline Surface::Halfedge Surface::Halfedge::prev() const {
        auto f = facet();
        int lh = id - f.halfedge();
        int n = f.size();
        return { m, f.halfedge((lh - 1 + n) % n) };
    }

    inline Surface::Halfedge Surface::Halfedge::opposite() const {
        assert(m.connected());
        assert(active());
        Halfedge result = { m, -1 }; // not found
        Halfedge cir = from().halfedge();
        while (true) {
            Halfedge candidate = cir.prev();
            if (cir.active()) {
                if (candidate.from() == to() && candidate.to() == from()) {
                    if (result == -1) result = candidate;
                    else return { m, -1 }; // found more than one
                }
                if (cir != *this && to() == cir.to())
                    return { m, -1 }; // the edge is non manifold
            }
            do {
                cir = m.conn->c2c[cir];
                if (cir == -1) return result;
            } while (!cir.active());
        }
    }

    inline Surface::Vertex Surface::Halfedge::from() const {
        auto f = facet();
        int lh = id - f.halfedge();
        return { m, m.vert(f, lh) };
    }

    inline Surface::Vertex Surface::Halfedge::to() const {
        auto f = facet();
        int lh = id - f.halfedge();
        int n = f.size();
        return { m, m.vert(f, (lh+1)%n) };
    }

    inline int Surface::Halfedge::id_in_facet() const {
        return id - m.corner(facet(), 0);
    }

    inline Surface::Halfedge::operator vec3() const {
        return {(vec3&)to() - (vec3&)from()};
    }

    inline Surface::Halfedge::operator Segment3() const {
        return {from(), to()};
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void Surface::Facet::deactivate() {
        assert(m.connected());
        m.conn->active[id] = false;
    }

    inline Surface::Vertex Surface::Facet::vertex(int lv) const {
        return { m, m.vert(id, lv) };
    }

    inline Surface::Halfedge Surface::Facet::halfedge(int lh) const {
        return { m, m.corner(id, lh) };
    }


    inline int Surface::Facet::size() const {
        return m.facet_size(id);
    }

    inline bool Surface::Facet::active() const {
        return id>=0 && (!m.connected() || m.conn->active[id]);
    }

    inline Surface::Facet::operator Triangle3() const {
        um_assert(size()==3);
        return { m.points[m.vert(id, 0)], m.points[m.vert(id, 1)], m.points[m.vert(id, 2)] };
    }

    inline Surface::Facet::operator Quad3() const {
        um_assert(size()==4);
        return { m.points[m.vert(id, 0)], m.points[m.vert(id, 1)], m.points[m.vert(id, 2)], m.points[m.vert(id, 3)] };
    }

    inline Surface::Facet::operator Poly3() const {
        std::vector<vec3> pts(size());
        for (int i = 0; i < size(); i++)
            pts[i] = m.points[m.vert(id, i)];
        return {pts};
    }

}

