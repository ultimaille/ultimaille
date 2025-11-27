#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <initializer_list>
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "attributes.h"
#include "pointset.h"
#include "primitive_geometry.h"
#include "syntactic-sugar/assert.h"

namespace UM {

    struct Surface { // polygonal mesh interface
        PointSet points{};
        std::vector<int> facets{};
        std::vector<std::weak_ptr<ContainerBase> > attr_facets{};
        std::vector<std::weak_ptr<ContainerBase> > attr_corners{};

////////////////////////////////////////////////////
//       _                 _                _     //
//      | | _____      __ | | _____   _____| |    //
//      | |/ _ \ \ /\ / / | |/ _ \ \ / / _ \ |    //
//      | | (_) \ V  V /  | |  __/\ V /  __/ |    //
//      |_|\___/ \_/\_/   |_|\___| \_/ \___|_|    //
//                                                //
////////////////////////////////////////////////////

        void delete_vertices(const std::vector<bool>& to_kill);
        virtual void delete_facets(const std::vector<bool>& to_kill);
        void delete_isolated_vertices();
        void resize_attrs();
        void compress_attrs(const std::vector<bool>& facets_to_kill);

        int nverts() const;
        int ncorners() const;

        virtual int nfacets() const = 0;
        virtual int facet_size(const int fi) const = 0;
        virtual int corner(const int fi, const int ci) const = 0;
        virtual int   vert(const int fi, const int lv) const = 0;
        virtual int&  vert(const int fi, const int lv) = 0;

        Surface() = default;
        Surface(const Surface& m) = delete;
        Surface(Surface&& m) = delete;
        Surface& operator=(const Surface& m) = delete;

//////////////////////////////////////////////////////////////////////
//                                      _   _       _ _             //
//       ___ ___  _ __  _ __   ___  ___| |_(_)_   _(_) |_ _   _     //
//      / __/ _ \| '_ \| '_ \ / _ \/ __| __| \ \ / / | __| | | |    //
//     | (_| (_) | | | | | | |  __/ (__| |_| |\ V /| | |_| |_| |    //
//      \___\___/|_| |_|_| |_|\___|\___|\__|_| \_/ |_|\__|\__, |    //
//                                                        |___/     //
//////////////////////////////////////////////////////////////////////

        struct Vertex;
        struct Halfedge;
        struct Facet;

        Vertex     vertex(int id) { return   Vertex(*this, id); }
        Halfedge halfedge(int id) { return Halfedge(*this, id); }
        Facet       facet(int id) { return    Facet(*this, id); }

        struct Connectivity {
            Surface& m;
            PointAttribute<int>  v2c;    // vertex to corner map
            CornerAttribute<int> c2f;    // corner to facet map
            CornerAttribute<int> c2c;    // corner to corner (sharing the same vertex) map. Consecutive maps form an unordered linked list terminated by -1.
            FacetAttribute<bool> active; // facets to keep after compacting

            Connectivity(Surface& m);
            void init();
            Surface::Facet create_facet(std::initializer_list<int> verts);
            Surface::Facet create_facet(int* verts, int size);

            void change_from(Surface::Halfedge he, int new_vertex_id); // TODO move it to the toolbox
        };

        std::unique_ptr<Connectivity> conn = {};
        inline bool connected() const { return conn != nullptr; }

        void connect();
        void disconnect();
        void compact(bool delete_isolated_vertices = true);

//////////////////////////////////////////////////////////////
//                  _           _ _   _                     //
//       _ __  _ __(_)_ __ ___ (_) |_(_)_   _____  ___      //
//      | '_ \| '__| | '_ ` _ \| | __| \ \ / / _ \/ __|     //
//      | |_) | |  | | | | | | | | |_| |\ V /  __/\__ \     //
//      | .__/|_|  |_|_| |_| |_|_|\__|_| \_/ \___||___/     //
//      |_|                                                 //
//////////////////////////////////////////////////////////////

        struct Primitive {
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

        struct Vertex : Primitive {
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

        struct Halfedge : Primitive {
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

            [[deprecated]] inline Segment3 geom();
            inline operator vec3() const;
            inline operator Segment3() const;

            auto iter_sector_halfedges() const;
        };

        struct Facet : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;

            void deactivate();

            Vertex vertex(int lv) const;
            Halfedge halfedge(int lh = 0) const;
            int size() const;
            bool active() const;

            template<typename T> [[deprecated]] T geom();
            operator Triangle3() const;
            operator Quad3() const;
            operator Poly3() const;

            auto iter_halfedges() const;
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

        int nfacets() const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int& vert(const int fi, const int lv);

        Triangles() = default;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Quads : Surface { // quad mesh implementation
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int& vert(const int fi, const int lv);

        Quads() = default;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Polygons : Surface { // polygonal mesh implementation
        std::vector<int> offset = { 0 };
        Polygons() = default;

        int create_facets(const int n, const int size);
        void delete_facets(const std::vector<bool>& to_kill);

        int nfacets()  const;
        int facet_size(const int fi) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int& vert(const int fi, const int lv);
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
        assert(0 == facets.size() % 3);
        return facets.size() / 3;
    }

    inline int Triangles::facet_size(const int) const {
        return 3;
    }

    inline int Triangles::corner(const int fi, const int ci) const {
        assert(ci >= 0 && ci < 3 && fi >= 0 && fi < nfacets());
        return fi * 3 + ci;
    }

    inline int Triangles::vert(const int fi, const int lv) const {
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 3);
        return facets[fi * 3 + lv];
    }

    inline int& Triangles::vert(const int fi, const int lv) {
        //      assert(conn==nullptr);
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 3);
        return facets[fi * 3 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Quads::nfacets() const {
        assert(0 == facets.size() % 4);
        return facets.size() / 4;
    }

    inline int Quads::facet_size(const int) const {
        return 4;
    }

    inline int Quads::corner(const int fi, const int ci) const {
        assert(ci >= 0 && ci < 4 && fi >= 0 && fi < nfacets());
        return fi * 4 + ci;
    }

    inline int Quads::vert(const int fi, const int lv) const {
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 4);
        return facets[fi * 4 + lv];
    }

    inline int& Quads::vert(const int fi, const int lv) {
        //      assert(conn==nullptr);
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 4);
        return facets[fi * 4 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Polygons::nfacets() const {
        return static_cast<int>(offset.size()) - 1;
    }

    inline int Polygons::facet_size(const int fi) const {
        assert(fi >= 0 && fi < nfacets());
        return offset[fi + 1] - offset[fi];
    }

    inline int Polygons::corner(const int fi, const int ci) const {
        assert(fi >= 0 && fi < nfacets());
        return offset[fi] + ci;
    }

    inline int Polygons::vert(const int fi, const int lv) const {
        assert(fi >= 0 && fi < nfacets());
        assert(lv >= 0 && lv < facet_size(fi));
        return facets[offset[fi] + lv];
    }

    inline int& Polygons::vert(const int fi, const int lv) {
        //      assert(conn==nullptr);
        assert(fi >= 0 && fi < nfacets());
        assert(lv >= 0 && lv < facet_size(fi));
        return facets[offset[fi] + lv];
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

    inline Segment3 Surface::Halfedge::geom() {
        return {from().pos(), to().pos()};
    }

    inline Surface::Halfedge::operator vec3() const {
        return {(vec3&)to() - (vec3&)from()};
    }

    inline Surface::Halfedge::operator Segment3() const {
        return {from().pos(), to().pos()};
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

    template<> inline Triangle3 Surface::Facet::geom() {
        um_assert(size()==3);
        return Triangle3(vertex(0).pos(), vertex(1).pos(), vertex(2).pos());
    }

    template<> inline Quad3 Surface::Facet::geom() {
        um_assert(size()==4);
        return Quad3(vertex(0).pos(), vertex(1).pos(), vertex(2).pos(), vertex(3).pos());
    }

    template<> inline Poly3 Surface::Facet::geom() {
        std::vector<vec3> pts(size());
        for (int i = 0; i < size(); i++)
            pts[i] = vertex(i).pos();

        return Poly3{pts};
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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

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

#endif //__SURFACE_H__
