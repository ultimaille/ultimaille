#ifndef __VOLUME_H__
#define __VOLUME_H__
#include <vector>
#include <array>
#include <memory>

#include "syntactic-sugar/assert.h"
#include "algebra/vec.h"
#include "attributes.h"
#include "pointset.h"
#include "volume_reference.h"
#include "volume_connectivity.h"
#include "primitive_geometry.h"

namespace UM {
    struct Volume;
    struct HalfEdgeHelper { // half-edge-like connectivity interface
        constexpr HalfEdgeHelper(const Volume &mesh) : m(mesh) {}

        constexpr int halfedge(const int cell, const int cell_facet, const int facet_he) const;
        int halfedge_from_verts(const int c, const int org, const int dst) const;

        int nhalfedges() const;
        constexpr int nhalfedges_per_cell() const;

        constexpr int           cell(const int he) const;
        constexpr int          facet(const int he) const;
        constexpr int         corner(const int he) const;
        constexpr int     cell_facet(const int he) const;
        constexpr int  cell_halfedge(const int he) const;
        constexpr int facet_halfedge(const int he) const;
        constexpr int           prev(const int he) const;
        constexpr int           next(const int he) const;
        constexpr int     opposite_f(const int he) const;

        vec3          geom(const int he) const;
        int           from(const int he) const;
        int             to(const int he) const;
        int opposite_c(const OppositeFacet &adj, const int he) const;
        const Volume &m;
    };

    // struct [[deprecated]] halfedge_around_edge_iter {
    //     const OppositeFacet  &of;
    //     const HalfEdgeHelper &heh;
    //     int start  = -1;
    //     int finish = -1;

    //     halfedge_around_edge_iter(const OppositeFacet &of, const int he);

    //     struct iterator {
    //         const OppositeFacet  &of;
    //         const HalfEdgeHelper &heh;
    //         int value;
    //         bool circ;

    //         void operator++() {
    //             value = of.opposite_c(heh.opposite_f(value));
    //             circ = false;
    //         }

    //         int operator*() const {
    //             return value;
    //         }

    //         bool operator!=(const iterator& rhs) const {
    //             return circ || value != rhs.value;
    //         }
    //     };

    //     iterator begin() const { return {of, heh, start, start==finish}; }
    //     iterator end()   const { return {of, heh, finish, false}; }
    // };

    struct Volume {
        enum CELL_TYPE { TETRAHEDRON=0, HEXAHEDRON=1, WEDGE=2, PYRAMID=3 };
        CELL_TYPE cell_type;

        // [[deprecated]] HalfEdgeHelper heh;
        PointSet points{};
        std::vector<int> cells{};

        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_cells{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_corners{};

        int  create_cells(const int n);
        void delete_cells(const std::vector<bool> &to_kill);
        void delete_vertices(const std::vector<bool> &to_kill);
        void delete_isolated_vertices();

        void resize_attrs();
        void compress_attrs(const std::vector<bool> &cells_to_kill);

        int nverts()   const;
        int ncells()   const;
        int nfacets()  const;
        int ncorners() const;
        constexpr int cell_from_facet (const int f) const;
        constexpr int cell_from_corner(const int c) const;
        int  vert(const int c, const int lv) const;
        int &vert(const int c, const int lv);

        constexpr int  nverts_per_cell() const;
        constexpr int nfacets_per_cell() const;
        constexpr int facet_size(const int f) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        constexpr int  facet(const int c, const int lf) const;
        constexpr int corner(const int c, const int lc) const;

        void clear() {
            points = {};
            cells  = {};
            attr_cells   = {};
            attr_facets  = {};
            attr_corners = {};
        }

        Volume(CELL_TYPE cell_type) : cell_type(cell_type) {}
        Volume(const Volume& m) { // TODO re-think copying policy
            um_assert(!m.points.size() && !m.cells.size());
        }
        Volume& operator=(const Volume& m) {
            clear();
            um_assert(!m.points.size() && !m.cells.size());
            return *this;
        }

        struct Connectivity {
            Volume &m;
            OppositeFacet oppf;
            HalfEdgeHelper heh;

            Connectivity(Volume &m);
            void reset();
        };

        std::unique_ptr<Connectivity> conn = {};
        inline bool connected() const { return conn != nullptr; }

        void connect();
        void disconnect();

        // TODO careful assert policy, esp. for the iterators

        struct Primitive {
            Primitive(Volume& m, int id) : m(m), id(id) {}
            Primitive(Primitive& p)  = default;
            Primitive(Primitive&& p) = default;

            Primitive& operator=(Primitive& p);
            Primitive& operator=(int i);

            operator int() const;
            operator int& ();
            bool active() const;

        protected:
            friend struct Volume;
            Volume& m;
            int id;
        };

        struct Vertex;
        struct Corner;
        struct Halfedge;
        struct Facet;
        struct Cell;

        struct Vertex :  Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Vertex(Vertex& v)  = default;
            Vertex(Vertex&& v) = default;
            Vertex& operator=(Vertex& v);

            vec3  pos() const;
            vec3& pos();
        };

        struct Corner : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Corner(Corner& v)  = default;
            Corner(Corner&& v) = default;
            Corner& operator=(Corner& v);

            int id_in_cell();
            Vertex vertex();
            Cell cell();
        };

        struct Halfedge : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Halfedge(Halfedge& he)  = default;
            Halfedge(Halfedge&& he) = default;
            Halfedge& operator=(Halfedge& he);

            Facet facet();
            int id_in_facet();
            Vertex from();
            Vertex to();
            Corner from_corner();
            Corner to_corner();
            Halfedge next();
            Halfedge prev();
            Cell cell();
            int id_in_cell();
            Halfedge opposite_f();
            Halfedge opposite_c();

            inline Segment3 geom();

            auto iter_CCW_around_edge();
        };

        struct Facet : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Facet(Facet& he)  = default;
            Facet(Facet&& he) = default;
            Facet& operator=(Facet& he);

            int nverts()     const;
            int ncorners()   const;
            int nhalfedges() const;

            Halfedge halfedge(int lh);
            //[[deprecated]]
            Vertex __vertex(int lv);
            Vertex vertex(int lv);
            Corner corner(int lc);

            Facet opposite();
            bool on_boundary();

            Cell cell();
            int id_in_cell();

            template<typename T> T geom();

            auto iter_halfedges();
        };

        struct Cell : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;
            Cell(Cell& he)  = default;
            Cell(Cell&& he) = default;
            Cell& operator=(Cell& he);

            int nfacets()    const;
            int ncorners()   const;
            int nverts()     const;
            int nhalfedges() const;

            Facet facet(int lf);
            Corner corner(int lc);
            Vertex vertex(int lv);
            Halfedge halfedge(int lh);

            template<typename T> T geom();

            auto iter_facets();
        };

        auto iter_vertices();
        auto iter_corners();
        auto iter_halfedges();
        auto iter_facets();
        auto iter_cells();
    };

   /*
    * EdgeGraph represents edges of a volumetric mesh
    *    -> Concept: An edge corresponds to a set of halfedges with same from/to vertices
    *    -> Design choice: Edges are manifold, so a non manifold edge is represented edges with same from/to vertices
    *   -> Benefits: A) use EdgeAttributes and B) gives acces to vertex neigborhood i.e. all primitives can reach its neigborgs
    */

    struct EdgeGraph : public PolyLine {
        EdgeGraph(Volume& m);
        PolyLine::Edge edge_from_halfedge(Volume::Halfedge h);
        Volume::Halfedge halfedge_from_edge(Edge e);

        std::vector<int> m_halfedge_from_edge;
        Volume& m;
    };


    inline Volume::Halfedge EdgeGraph::halfedge_from_edge(Edge e) { return { m, m_halfedge_from_edge[e] }; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Tetrahedra : Volume {
        Tetrahedra() : Volume(Volume::TETRAHEDRON) {}
    };

    struct Hexahedra : Volume {
        Hexahedra() : Volume(Volume::HEXAHEDRON) {}
    };

    struct Wedges : Volume {
        Wedges() : Volume(Volume::WEDGE) {}
    };

    struct Pyramids : Volume {
        Pyramids() : Volume(Volume::PYRAMID) {}
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // these implementations are here and not in the .cpp because all inline functions must be available in all translation units //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Volume::nverts() const {
        return points.size();
    }

    inline int Volume::ncorners() const {
        return cells.size();
    }

    inline int Volume::ncells() const {
        return cells.size()/nverts_per_cell();
    }

    inline int Volume::nfacets() const {
        return ncells()*nfacets_per_cell();
    }

    inline constexpr int Volume::cell_from_facet(const int f) const {
        assert(f>=0 && f<nfacets());
        return f/nfacets_per_cell();
    }

    inline constexpr int Volume::cell_from_corner(const int c) const {
        assert(c>=0 && c<ncorners());
        return c/nverts_per_cell();
    }

    inline int Volume::vert(const int c, const int lv) const {
        assert(c>=0 && c<ncells() && lv>=0 && lv<nverts_per_cell());
        return cells[corner(c, lv)];
    }

    inline int &Volume::vert(const int c, const int lv) {
        assert(c>=0 && c<ncells() && lv>=0 && lv<nverts_per_cell());
        return cells[corner(c, lv)];
    }

    inline constexpr int Volume::nverts_per_cell() const {
        return reference_cells[cell_type].nverts();
    }

    inline constexpr int Volume::nfacets_per_cell() const {
        return reference_cells[cell_type].nfacets();
    }

    inline constexpr int Volume::facet_size(const int f) const {
        assert(f>=0 && f<nfacets());
        return reference_cells[cell_type].facet_size(f % nfacets_per_cell());
    }

    inline int Volume::facet_vert(const int c, const int lf, const int lv) const {
        // TODO check assert ! facet_size(lf) but facet_size expect a facet index, not local index !
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell() && lv>=0 && lv<facet_size(lf));
        return vert(c, reference_cells[cell_type].vert(lf, lv));
    }

    inline constexpr int Volume::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        return c*nfacets_per_cell() + lf;
    }

    inline constexpr int Volume::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<nverts_per_cell());
        return c*nverts_per_cell() + lc;
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

    inline constexpr int HalfEdgeHelper::nhalfedges_per_cell() const {
        return reference_cells[m.cell_type].ncorners();
    }

    // global halfedge id from local id
    inline constexpr int HalfEdgeHelper::halfedge(const int cell, const int cell_facet, const int facet_he) const {
//      assert(cell>=0 && cell<m.ncells());
        assert(cell_facet>=0 && cell_facet<m.nfacets_per_cell());
        assert(facet_he>=0 && facet_he<=m.facet_size(cell_facet));
        return cell*nhalfedges_per_cell() + reference_cells[m.cell_type].corner(cell_facet, facet_he);
    }

    // global cell id
    inline constexpr int HalfEdgeHelper::cell(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return he / nhalfedges_per_cell();
    }

    // global facet id
    inline constexpr int HalfEdgeHelper::facet(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return m.facet(cell(he), cell_facet(he));
    }

    // global corner id
    inline constexpr int HalfEdgeHelper::corner(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return cell(he) * m.nverts_per_cell() + reference_cells[m.cell_type].from(cell_halfedge(he));
    }

    // local facet id
    inline constexpr int HalfEdgeHelper::cell_facet(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return reference_cells[m.cell_type].facet(cell_halfedge(he));
    }

    // local halfedge id
    inline constexpr int HalfEdgeHelper::cell_halfedge(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return he % nhalfedges_per_cell();
    }

    // local halfedge id
    inline constexpr int HalfEdgeHelper::facet_halfedge(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return cell_halfedge(he) - reference_cells[m.cell_type].corner(cell_facet(he), 0);
    }

    // global cell halfedge id
    inline constexpr int HalfEdgeHelper::prev(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        if (facet_halfedge(he)>0) return he - 1;
        const int size = m.facet_size(cell_facet(he));
        return he + size - 1;
    }

    // global cell halfedge id
    inline constexpr int HalfEdgeHelper::next(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        const int size = m.facet_size(cell_facet(he));
        if (facet_halfedge(he)<size-1) return he + 1;
        return he - size + 1;
    }

    inline constexpr int HalfEdgeHelper::opposite_f(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return nhalfedges_per_cell()*cell(he) + reference_cells[m.cell_type].opposite(cell_halfedge(he));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int HalfEdgeHelper::nhalfedges() const {
        return m.ncells() * nhalfedges_per_cell();
    }

    inline vec3 HalfEdgeHelper::geom(const int he) const {
        return m.points[to(he)] - m.points[from(he)];
    }

    // global vertex id
    inline int HalfEdgeHelper::from(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return m.facet_vert(cell(he), cell_facet(he), facet_halfedge(he));
    }

    // global vertex id
    inline int HalfEdgeHelper::to(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return from(next(he));
    }

    // inline int HalfEdgeHelper::opposite_c(const OppositeFacet &adj, const int he) const {
    //     return adj.opposite_c(he);
    // }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Volume::Primitive::operator int() const {
        return id;
    }

    inline Volume::Primitive::operator int& () {
        return id;
    }

    inline bool Volume::Primitive::active() const {
        return id>=0;
    }

    inline Volume::Primitive& Volume::Primitive::operator=(Volume::Primitive& p) {
        assert(&m == &p.m);
        id = p.id;
        return *this;
    }

    inline Volume::Primitive& Volume::Primitive::operator=(int i) {
        id = i;
        return *this;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Volume::Vertex& Volume::Vertex::operator=(Volume::Vertex& v) {
        Primitive::operator=(v);
        return *this;
    }

    inline vec3  Volume::Vertex::pos() const {
        return m.points[id];
    }

    inline vec3& Volume::Vertex::pos() {
        return m.points[id];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Volume::Corner& Volume::Corner::operator=(Volume::Corner& v) {
        Primitive::operator=(v);
        return *this;
    }

    inline Volume::Vertex Volume::Corner::vertex() {
        return { m, m.vert(id / m.nverts_per_cell(), id % m.nverts_per_cell()) };
    }

    inline Volume::Cell Volume::Corner::cell() {
         return { m, id / m.nverts_per_cell() };
    }

    inline int Volume::Corner::id_in_cell() {
        return id % m.nverts_per_cell();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Volume::Halfedge& Volume::Halfedge::operator=(Volume::Halfedge& v) {
        Primitive::operator=(v);
        return *this;
    }

    inline Volume::Vertex Volume::Halfedge::from() {
        assert(m.connected());
        return { m, m.conn->heh.from(id) };
    }

    inline Volume::Vertex Volume::Halfedge::to() {
        assert(m.connected());
        return { m, m.conn->heh.to(id) };
    }

    inline Volume::Corner Volume::Halfedge::from_corner() {
        assert(m.connected());
        return { m, m.conn->heh.corner(id) };
    }

    inline Volume::Corner Volume::Halfedge::to_corner() {
        assert(m.connected());
        return { m, next().from_corner() };
    }

    inline Volume::Halfedge Volume::Halfedge::next() {
        assert(m.connected());
        return { m, m.conn->heh.next(id) };
    }

    inline Volume::Halfedge Volume::Halfedge::prev() {
        assert(m.connected());
        return { m, m.conn->heh.prev(id) };
    }

    inline Volume::Facet Volume::Halfedge::facet() {
        assert(m.connected());
        return { m, m.conn->heh.facet(id) };
    }

    inline int Volume::Halfedge::id_in_facet() {
        assert(m.connected());
        return id - facet().halfedge(0);
    }

    inline Volume::Cell Volume::Halfedge::cell() {
        assert(m.connected());
        return { m, m.conn->heh.cell(id) };
    }

    inline int Volume::Halfedge::id_in_cell() {
        assert(m.connected());
        return id - m.conn->heh.cell(id)*m.conn->heh.nhalfedges_per_cell();
    }

    inline Volume::Halfedge Volume::Halfedge::opposite_f() {
        assert(m.connected());
        return { m, m.conn->heh.opposite_f(id) };
    }

    inline Volume::Halfedge Volume::Halfedge::opposite_c() {
        assert(m.connected());
        Facet oppf = facet().opposite();
        if (!oppf.active()) return { m, -1 };
        for (int lv=0; lv<oppf.nhalfedges(); lv++) {
            Halfedge res = oppf.halfedge(lv);
            if (res.from() == to() && res.to() == from()) return res;
        }
        um_assert(false);
        return { m, -1 };
    }

    inline Segment3 Volume::Halfedge::geom() {
        return {from().pos(), to().pos()};
    }

    inline auto Volume::Halfedge::iter_CCW_around_edge() {
        struct iterator {
            Halfedge ref;
            Halfedge data;
            void operator++() {
                data = data.opposite_f().opposite_c();
                if (data == ref) data.id = -1;
            }
            bool operator!=(iterator& rhs) { return data != rhs.data; }
            Halfedge& operator*() { return data; }
        };
        struct wrapper {
            Halfedge ref;
            auto begin() {
                if (!ref.active())return iterator{ ref,ref };
                Halfedge org = ref;
                do {
                    Halfedge opp = ref.opposite_c();
                    if (!opp.active())
                        return iterator{ ref , ref };
                    ref = opp.opposite_f();
                } while (org != ref);
                return iterator{ ref , ref };
            }
            auto end() { return iterator{ ref , Halfedge(ref.m,-1) }; }
        };
        return wrapper{ *this };
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Volume::Facet& Volume::Facet::operator=(Volume::Facet& v) {
        Primitive::operator=(v);
        return *this;
    }

    inline int Volume::Facet::nverts() const {
        return m.facet_size(id);
    }

    inline int Volume::Facet::ncorners() const {
        return m.facet_size(id);
    }

    inline int Volume::Facet::nhalfedges() const {
        return m.facet_size(id);
    }

    inline Volume::Halfedge Volume::Facet::halfedge(int i) {
        assert(m.connected());
        return { m, m.conn->heh.halfedge(cell(), id_in_cell(), i) };
    }

    // No need to be connected anymore
    inline Volume::Vertex Volume::Facet::__vertex(int i) {
        assert(m.connected());
        return { m, halfedge(i).from() };
    }

    inline Volume::Vertex Volume::Facet::vertex(int i) {
        return Volume::Vertex(m, m.facet_vert(cell(), id_in_cell(), i));
    }

            

    inline Volume::Corner Volume::Facet::corner(int i) {
        assert(m.connected());
        return { m, halfedge(i).from_corner() };
    }

    inline Volume::Facet Volume::Facet::opposite() {
        assert(m.connected());
        return { m, m.conn->oppf.adjacent[id] };
    }

    inline bool Volume::Facet::on_boundary() {
        return !opposite().active();
    }

    inline Volume::Cell Volume::Facet::cell() {
        //assert(m.connected());
        return { m, m.cell_from_facet(id) };
    }

    inline int Volume::Facet::id_in_cell() {
        return id % m.nfacets_per_cell();
    }

    inline auto Volume::Facet::iter_halfedges() {
        struct iterator {
            Halfedge data;
            void operator++() { int f = data.facet(); ++(data.id); if (data.id % data.m.facet_size(f) == 0) data.id = -1; }
            bool operator!=(iterator& rhs) { return data != rhs.data; }
            Halfedge& operator*() { return data; }
        };
        struct wrapper {
            Facet f;
            auto begin() { return iterator{ f.halfedge(0) }; }
            auto end()   { return iterator{ Halfedge(f.m,-1) }; }
        };
        return wrapper{ *this };
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Volume::Cell& Volume::Cell::operator=(Volume::Cell& v) {
        Primitive::operator=(v);
        return *this;
    }

    inline int Volume::Cell::nverts() const {
        return m.nverts_per_cell();
    }

    inline int Volume::Cell::nfacets() const {
        return m.nfacets_per_cell();
    }

    inline int Volume::Cell::ncorners() const {
        return m.nverts_per_cell();
    }

    inline int Volume::Cell::nhalfedges() const {
        assert(m.connected());
        return m.conn->heh.nhalfedges_per_cell();
    }

    inline Volume::Vertex Volume::Cell::vertex(int lv) {
        //assert(m.connected());
        return { m, m.vert(id, lv) };
    }

    inline Volume::Corner Volume::Cell::corner(int lc) {
//        assert(m.connected());
        return { m, m.corner(id, lc) };
    }

    inline Volume::Halfedge Volume::Cell::halfedge(int lh) {
        assert(m.connected());
        return { m, m.conn->heh.nhalfedges_per_cell()*id + lh };
    }

    inline Volume::Facet Volume::Cell::facet(int lf) {
//        assert(m.connected());
        return { m, m.nfacets_per_cell()*id + lf };
    }

    inline auto Volume::Cell::iter_facets() {
        struct iterator {
            Facet data;
            void operator++() { ++(data.id); }
            bool operator!=(iterator& rhs) { return data != rhs.data; }
            Facet& operator*() { return data; }
        };
        struct wrapper {
            Cell c;
            auto begin() { return iterator{ c.facet(0) }; }
            auto end()   { return iterator{ c.facet(c.nfacets()) }; }
        };
        return wrapper{ *this };
    }

    template<> inline Tetrahedron Volume::Cell::geom() {
        um_assert(nfacets()==4 && nverts()==4);
        return Tetrahedron(vertex(0).pos(), vertex(1).pos(), vertex(2).pos(), vertex(3).pos());
    }

    template<> inline Pyramid Volume::Cell::geom() {
        um_assert(nfacets()==5 && nverts()==5);
        return Pyramid(vertex(0).pos(), vertex(1).pos(), vertex(2).pos(), vertex(3).pos(), vertex(4).pos());
    }

    template<> inline Hexahedron Volume::Cell::geom() {
        um_assert(nfacets()==6 && nverts()==8);
        return Hexahedron(vertex(0).pos(), vertex(1).pos(), vertex(2).pos(), vertex(3).pos(), vertex(4).pos(), vertex(5).pos(), vertex(6).pos(), vertex(7).pos());
    }

    template<> inline Wedge Volume::Cell::geom() {
        um_assert(nfacets()==5 && nverts()==6);
        return Wedge(vertex(0).pos(), vertex(1).pos(), vertex(2).pos(), vertex(3).pos(), vertex(4).pos(), vertex(5).pos());
    }

    template<> inline Polyhedron Volume::Cell::geom() {
        std::vector<vec3> pts(nverts());
        for (int i = 0; i < nverts(); i++)
            pts[i] = vertex(i).pos();

        return Polyhedron{pts};
    }

    template<> inline Triangle3 Volume::Facet::geom() {
        um_assert(nverts()==3);
        return Triangle3(vertex(0).pos(), vertex(1).pos(), vertex(2).pos());
    }

    template<> inline Quad3 Volume::Facet::geom() {
        um_assert(nverts()==4);
        return Quad3(vertex(0).pos(), vertex(1).pos(), vertex(2).pos(), vertex(3).pos());
    }

    template<> inline Poly3 Volume::Facet::geom() {
        // TODO replace m.facet_size(id) => size() but we have to add size() on volume facet
        int size = m.facet_size(id);
        std::vector<vec3> pts(size);
        for (int i = 0; i < size; i++)
            pts[i] = vertex(i).pos();

        return Poly3{pts};
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline auto Volume::iter_vertices() {
        struct iterator {
            Vertex v;
            void operator++() { ++v.id; }
            bool operator!=(iterator& rhs) { return v != rhs.v; }
            Vertex& operator*() { return v; }
        };
        struct wrapper {
            Volume& m;
            auto begin() { return iterator{ Vertex(m,0) }; }
            auto end() { return iterator{ Vertex(m,m.nverts()) }; }
        };
        return wrapper{ *this };
    }

    inline auto Volume::iter_corners() {
        struct iterator {
            Corner c;
            void operator++() { ++c.id; }
            bool operator!=(iterator& rhs) { return c != rhs.c; }
            Corner& operator*() { return c; }
        };
        struct wrapper {
            Volume& m;
            auto begin() { return iterator{ Corner(m,0) }; }
            auto end() { return iterator{ Corner(m,m.ncorners()) }; }
        };
        return wrapper{ *this };
    }

    inline auto Volume::iter_halfedges() {
        assert(connected());
        struct iterator {
            Halfedge data;
            void operator++() { ++data.id; }
            bool operator!=( iterator& rhs)  { return data != rhs.data; }
            Halfedge& operator*()  { return data; }
        };
        struct wrapper {
            Volume& m;
            auto begin() { return iterator{ Halfedge(m,0) }; }
            auto end() { return iterator{ Halfedge(m,m.conn->heh.nhalfedges()) }; }
        };
        return wrapper{ *this };
    }

    inline auto Volume::iter_facets() {
        struct iterator {
            Facet data;
            void operator++() { ++(data.id); }
            bool operator!=( iterator& rhs)  { return data != rhs.data; }
             Facet& operator*()  { return data; }
        };
        struct wrapper {
            Volume& m;
            auto begin() { return iterator{ Facet(m,0) }; }
            auto end() { return iterator{ Facet(m,m.nfacets()) }; }
        };
        return wrapper{ *this };
    }

    inline auto Volume::iter_cells() {
        struct iterator {
            Cell data;
            void operator++() { ++data.id; }
            bool operator!=( iterator& rhs)  { return data != rhs.data; }
            Cell& operator*()  { return data; }
        };
        struct wrapper {
            Volume& m;
            auto begin() { return iterator{ Cell(m,0) }; }
            auto end() { return iterator{ Cell(m,m.ncells()) }; }
        };
        return wrapper{ *this };
    }

}

#endif //__VOLUME_H__
