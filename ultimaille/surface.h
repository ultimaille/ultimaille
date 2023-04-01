#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "attributes.h"
#include "pointset.h"
#include "syntactic-sugar/assert.h"

namespace UM {
    struct Surface { // polygonal mesh interface
        PointSet points{};
        std::vector<int> facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_corners{};

        void delete_vertices(const std::vector<bool> &to_kill);
        virtual void delete_facets(const std::vector<bool> &to_kill);
        void delete_isolated_vertices();
        void resize_attrs();
        void compress_attrs(const std::vector<bool> &facets_to_kill);

        int nverts() const;
        int ncorners() const;

        virtual int nfacets() const = 0;
        virtual int facet_size(const int fi) const = 0;
        virtual int corner(const int fi, const int ci) const = 0;
        virtual int   vert(const int fi, const int lv) const = 0;
        virtual int  &vert(const int fi, const int lv)       = 0;

        virtual void clear() {
            points = {};
            attr_facets  = {};
            attr_corners = {};
            disconnect();
        }

        Surface() = default;
        Surface(const Surface& m) {
            um_assert(!m.points.size() && !m.facets.size());
        }
        Surface& operator=(const Surface& m) {
            clear();
            um_assert(!m.points.size() && !m.facets.size());
            return *this;
        }

        struct Util {
            const Surface &m;
            vec3 bary_verts(const int f) const;
        } util = {*this};

        struct Connectivity {
            Connectivity(Surface &mesh);
            PointAttribute<int>  v2c; // vertex to corner
            CornerAttribute<int> c2f; // corner to facet
            CornerAttribute<int> c2c; // corner to next corner sharing the same vertex (unordered)
            FacetAttribute<bool> active; // facets to keep after compacting
        };

        std::unique_ptr<Connectivity> conn = {};
        inline bool connected() const { return conn!=nullptr; }

        void connect();
        void disconnect();
        void compact(bool delete_isolated_vertices = true);

        struct Primitive {
            Primitive(Surface &m, int id);
            Primitive(Primitive &p) = default;

            Primitive& operator=(Primitive &p);
            Primitive& operator=(int i);

            operator int  () const;
            operator int& ();

            protected:
            friend struct Surface;

            Surface &m;
            int id;
        };

        struct Halfedge;
        struct Vertex : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Vertex(Vertex &v) = default;
            Vertex& operator=(const Vertex& v);

            vec3  pos() const;
            vec3 &pos();
            Halfedge halfedge();
            bool on_boundary();

            auto iter_halfedges();
        };

        struct Facet;
        struct Halfedge : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Halfedge(Halfedge &v) = default;
            Halfedge& operator=(const Halfedge& he);

            bool active();
            Facet facet();

            Halfedge next();
            Halfedge prev();
            Halfedge opposite();

            friend struct Vertex;
            Vertex from();
            Vertex to();

            auto iter_sector_halfedges();
        };

        struct Facet : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Facet(Facet &v) = default;
            Facet& operator=(const Facet& f);

            Vertex vertex(int lv);
            Halfedge halfedge(int lh = 0);
            int size();
            bool active();

            auto iter_halfedges();
        };

        auto iter_vertices();
        auto iter_halfedges();
        auto iter_facets();
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // these implementations are here and not in the .cpp because all inline functions must be available in all translation units //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Triangles : Surface { // simplicial mesh implementation
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        Triangles() = default;
        Triangles(const Triangles& m) = default;
        Triangles& operator=(const Triangles& m) {
            Surface::operator=(m);
            return *this;
        }

        struct Util : Surface::Util {
            double unsigned_area(const int f) const;
            void project(const int t, vec2& z0, vec2& z1, vec2& z2) const;
            vec3 normal(const int f) const;
        } util = {*this};
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Quads : Surface { // quad mesh implementation
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        Quads() = default;
        Quads(const Quads& m) = default;
        Quads& operator=(const Quads& m) {
            Surface::operator=(m);
            return *this;
        }

        struct Util : Surface::Util {
            double unsigned_area(const int f) const;
            vec3 normal(const int f) const;
        } util = {*this};
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Polygons : Surface { // polygonal mesh implementation
        std::vector<int> offset = {0};
        Polygons() = default;
        Polygons(const Polygons& m) = default;
        Polygons& operator=(const Polygons& m) {
            clear();
            um_assert(!m.points.size() && !m.facets.size() && 1==m.offset.size());
            return *this;
        }

        int create_facets(const int n, const int size);
        void delete_facets(const std::vector<bool> &to_kill);

        virtual void clear() {
            Surface::clear();
            offset = {0};
        }

        int nfacets()  const;
        int facet_size(const int fi) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        struct Util : Surface::Util {
        } util = {*this};
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Surface::nverts() const {
        return points.size();
    }

    inline int Surface::ncorners() const {
        return facets.size();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Triangles::nfacets() const {
        assert(0==facets.size()%3);
        return facets.size()/3;
    }

    inline int Triangles::facet_size(const int) const {
        return 3;
    }

    inline int Triangles::corner(const int fi, const int ci) const {
        assert(ci>=0 && ci<3 && fi>=0 && fi<nfacets());
        return fi*3 + ci;
    }

    inline int Triangles::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
        return facets[fi*3 + lv];
    }

    inline int &Triangles::vert(const int fi, const int lv) {
//      assert(conn==nullptr);
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
        return facets[fi*3 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Quads::nfacets() const {
        assert(0==facets.size()%4);
        return facets.size()/4;
    }

    inline int Quads::facet_size(const int) const {
        return 4;
    }

    inline int Quads::corner(const int fi, const int ci) const {
        assert(ci>=0 && ci<4 && fi>=0 && fi<nfacets());
        return fi*4 + ci;
    }

    inline int Quads::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<4);
        return facets[fi*4 + lv];
    }

    inline int &Quads::vert(const int fi, const int lv) {
//      assert(conn==nullptr);
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<4);
        return facets[fi*4 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Polygons::nfacets() const {
        return static_cast<int>(offset.size())-1;
    }

    inline int Polygons::facet_size(const int fi) const {
        assert(fi>=0 && fi<nfacets());
        return offset[fi+1]-offset[fi];
    }

    inline int Polygons::corner(const int fi, const int ci) const {
        assert(fi>=0 && fi<nfacets());
        return offset[fi]+ci;
    }

    inline int Polygons::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets());
        assert(lv>=0 && lv<facet_size(fi));
        return facets[offset[fi]+lv];
    }

    inline int &Polygons::vert(const int fi, const int lv) {
//      assert(conn==nullptr);
        assert(fi>=0 && fi<nfacets());
        assert(lv>=0 && lv<facet_size(fi));
        return facets[offset[fi]+lv];
    }

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

    inline Surface::Primitive::Primitive(Surface &m, int id) : m(m), id(id) {}

    inline Surface::Primitive& Surface::Primitive::operator=(Surface::Primitive &p) {
        assert(this == &p);
        id = p.id;
        return *this;
    }

    inline Surface::Primitive& Surface::Primitive::operator=(int i) {
        id = i;
        return *this;
    }

    inline Surface::Primitive::operator int  () const { return id; }
    inline Surface::Primitive::operator int& ()       { return id; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Surface::Vertex& Surface::Vertex::operator=(const Surface::Vertex& v) {
        Primitive::operator=(v);
        return *this;
    }

    inline vec3  Surface::Vertex::pos() const { return m.points[id]; }
    inline vec3 &Surface::Vertex::pos()       { return m.points[id]; }

    inline Surface::Halfedge Surface::Vertex::halfedge() {
        assert(m.connected());
        Surface::Halfedge res{m, m.conn->v2c[id]};
        while (res>=0 && !res.active()) {
            res = m.conn->c2c[res];
            if (res == m.conn->v2c[id])
                return {m, -1};
        }
        return res;
    }

    inline bool Surface::Vertex::on_boundary() {
        assert(m.connected());
        Surface::Halfedge cir = halfedge();
        if (cir<0) return false;
        do {
            if (cir.opposite() == -1) return true;
            cir = m.conn->c2c[cir];
        } while (cir != m.conn->v2c[id]);
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Surface::Halfedge& Surface::Halfedge::operator=(const Surface::Halfedge& he) {
        Primitive::operator=(he);
        return *this;
    }

    inline bool Surface::Halfedge::active() {
        return facet().active();
    }

    inline Surface::Facet Surface::Halfedge::facet() {
        assert(m.connected());
        return { m, m.conn->c2f[id] };
    }

    inline Surface::Halfedge Surface::Halfedge::next() {
        auto f = facet();
        int lh = id - f.halfedge();
        int n = f.size();
        return { m, f.halfedge((lh+1)%n) };
    }

    inline Surface::Halfedge Surface::Halfedge::prev() {
        auto f = facet();
        int lh = id - f.halfedge();
        int n = f.size();
        return { m, f.halfedge((lh-1+n)%n) };
    }

    inline Surface::Halfedge Surface::Halfedge::opposite() {
        assert(m.connected());
        assert(active());
        Halfedge cir = *this;
        Halfedge result = {m, -1}; // not found
        do {
            Halfedge candidate = cir.prev();
            if (cir.active()) {
                if (candidate.from() == to() && candidate.to() == from()) {
                    if (result == -1) result = candidate;
                    else return {m, -1}; // found more than one
                }
                if (cir != *this && to() == cir.to())
                    return {m, -1}; // the edge is non manifold
            }
            cir = cir.m.conn->c2c[cir];
        } while (cir != *this);
        return result;
    }

    inline Surface::Vertex Surface::Halfedge::from() {
        auto f = facet();
        int lh = id - f.halfedge();
        return { m, m.vert(f, lh) };
    }

    inline Surface::Vertex Surface::Halfedge::to() {
        auto f = facet();
        int lh = id - f.halfedge();
        int n = f.size();
        return { m, m.vert(f, (lh+1)%n) };
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Surface::Facet& Surface::Facet::operator=(const Surface::Facet& f) {
        Primitive::operator=(f);
        return *this;
    }

    inline Surface::Vertex Surface::Facet::vertex(int lv) {
        return { m, m.vert(id, lv) };
    }

    inline Surface::Halfedge Surface::Facet::halfedge(int lh) {
        return { m, m.corner(id, lh) };
    }

    inline int Surface::Facet::size() {
        return m.facet_size(id);
    }

    inline bool Surface::Facet::active() {
        return !m.connected() || m.conn->active[id];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO: is it really reasonable to inline these functions?

    inline auto Surface::iter_vertices() {
        struct iterator {
            Vertex v;
            void operator++() { ++v.id; }
            bool operator!=(const iterator& rhs) const { return v != rhs.v; }
            Vertex& operator*() { return v; }
        };
        struct wrapper {
            Surface& m;
            auto begin() { return iterator{ {m,0} }; }
            auto end() { return iterator{ {m,m.nverts()} }; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::iter_halfedges() {
        struct iterator {
            Halfedge h;
            void operator++() {
                ++h.id;
                if (h.m.connected())
                    while (h.id < h.m.ncorners() && !h.active()) ++h.id;
            }
            bool operator!=(const iterator& rhs) const { return h != rhs.h; }
            Halfedge& operator*() { return h; }
        };
        struct wrapper {
            Surface& m;
            auto begin() {
                iterator res{ {m,0} };
                if (m.connected())
                    if (m.ncorners() > 0 && !res.h.active()) ++res;
                return res;
            }
            auto end() { return iterator{ {m,m.ncorners()} }; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::iter_facets() {
        struct iterator {
            Facet f;
            void operator++() {
                ++f.id;
                if (f.m.connected())
                    while (f.id < f.m.nfacets() && !f.active()) ++f.id;
            }
            bool operator!=(const iterator& rhs) const { return f != rhs.f; }
            Facet& operator*() { return f; }
        };
        struct wrapper {
            Surface& m;
            auto begin() {
                iterator res{ {m,0} };
                if (m.connected())
                    if (m.nfacets() > 0 && !res.f.active()) ++res;
                return res;
            }
            auto end() { return iterator{ {m,m.nfacets()} }; }
        };
        return wrapper{ *this };
    }

    inline auto Surface::Vertex::iter_halfedges() {
        assert(m.connected());
        struct iterator {
            Surface::Halfedge he;
            void operator++() {
                const auto &c2c = he.m.conn->c2c;
                he = c2c[he];
                while (!he.active()) he = c2c[he];
                if (he == he.from().halfedge()) he = -1;
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

    inline auto Surface::Halfedge::iter_sector_halfedges() {
        struct iterator {
            Surface::Halfedge h, seed;
//          iterator(Surface::Halfedge seed) : h(seed), seed(seed) {}
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

    inline auto Surface::Facet::iter_halfedges() {
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

#endif //__SURFACE_H__

