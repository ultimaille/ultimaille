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
#include "primitive_geometry.h"
#include "polyline.h"

namespace UM {
    struct Volume;

    struct Volume {
        enum CELL_TYPE { TETRAHEDRON=0, HEXAHEDRON=1, WEDGE=2, PYRAMID=3 };
        CELL_TYPE cell_type;

        PointSet points{};
        std::vector<int> cells{};

        std::vector<std::weak_ptr<ContainerBase> > attr_cells{};
        std::vector<std::weak_ptr<ContainerBase> > attr_facets{};
        std::vector<std::weak_ptr<ContainerBase> > attr_corners{};

        int  create_cells(const int n);
        void delete_vertices(const std::vector<bool> &to_kill); // TODO invocable
        void delete_isolated_vertices();

        template <typename T> void delete_cells(const T &to_kill);

        void resize_attrs();

        int nverts()     const;
        int ncells()     const;
        int nfacets()    const;
        int ncorners()   const;
        int nhalfedges() const;

        struct Vertex; struct Corner; struct Halfedge; struct Facet; struct Cell;

        Vertex     vertex(int id) const { return   Vertex(*this, id); }
        Corner     corner(int id) const { return   Corner(*this, id); }
        Halfedge halfedge(int id) const { return Halfedge(*this, id); }
        Facet       facet(int id) const { return    Facet(*this, id); }
        Cell         cell(int id) const { return     Cell(*this, id); }

        constexpr int     nverts_per_cell() const;
        constexpr int    nfacets_per_cell() const;
        constexpr int nhalfedges_per_cell() const;
        constexpr int facet_size(const int f) const;

        [[deprecated]] constexpr int cell_from_facet (const int f) const;
        [[deprecated]] constexpr int cell_from_corner(const int c) const;
        [[deprecated]] constexpr int  facet(const int c, const int lf) const;
        [[deprecated]] constexpr int corner(const int c, const int lc) const;

        int facet_vert(const int c, const int lf, const int lv) const;
        int  vert(const int c, const int lv) const;
        int &vert(const int c, const int lv);

        Volume(CELL_TYPE cell_type) : cell_type(cell_type) {}
        Volume(const Volume& m)            = delete;
        Volume(Volume&& m)                 = delete;
        Volume& operator=(const Volume& m) = delete;

        struct Connectivity {
            Volume &m;
            CellFacetAttribute<int> adjacent;

            Connectivity(Volume &m);
            void reset();
        };

        std::unique_ptr<Connectivity> conn = {};
        inline bool connected() const { return conn != nullptr; }

        void connect();
        void disconnect();

        // TODO careful assert policy, esp. for the iterators

        struct Primitive {
            Primitive(const Volume& m, int id) : m(m), id(id) {}
            Primitive(Primitive& p)  = default;
            Primitive(Primitive&& p) = default;
            Primitive(const Primitive& p)  = default;

            Primitive& operator=(const Primitive&& p) noexcept;
            Primitive& operator=(Primitive& p);
            Primitive& operator=(int i);

            operator int() const;
            operator int& ();
            bool active() const;

        protected:
            friend struct Volume;
            const Volume& m;
            int id;
        };

        struct Vertex : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;

            [[deprecated]] vec3  pos() const;
            [[deprecated]] vec3& pos();

            inline operator const vec3&() const;
        };

        struct Corner : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;

            int id_in_cell()    const;
            Vertex vertex()     const;
            Cell cell()         const;
            Halfedge halfedge() const;

//          auto iter_halfedges() const; // TODO
        };

        struct Halfedge : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;

            int id_in_facet()     const;
            int id_in_cell()      const;
            Vertex from()         const;
            Vertex to()           const;
            Corner from_corner()  const;
            Corner to_corner()    const;
            Halfedge next()       const;
            Halfedge prev()       const;
            Halfedge opposite_f() const;
            Halfedge opposite_c() const;
            Facet facet()         const;
            Cell cell()           const;

            operator Segment3()   const;

            auto iter_CCW_around_edge();
        };

        struct Facet : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;

            [[deprecated]] int nverts()     const;
            [[deprecated]] int ncorners()   const;
            [[deprecated]] int nhalfedges() const;
            int size() const;

            bool on_boundary() const;

            int id_in_cell()          const;
            Vertex vertex(int lv)     const;
            Halfedge halfedge(int lh) const;
            Corner corner(int lc)     const;
            Facet opposite()          const;
            Cell cell()               const;

            operator Triangle3() const;
            operator Quad3()     const;
            operator Poly3()     const;

            auto iter_halfedges();
        };

        struct Cell : Primitive {
            using Primitive::Primitive;
            using Primitive::operator=;

            int nfacets()    const;
            int ncorners()   const;
            int nverts()     const;
            int nhalfedges() const;

            Facet facet(int lf)       const;
            Corner corner(int lc)     const;
            Vertex vertex(int lv)     const;
            Halfedge halfedge(int lh) const;

            operator Tetrahedron() const;
            operator Hexahedron()  const;
            operator Pyramid()     const;
            operator Wedge()       const;

            auto iter_halfedges();
            auto iter_facets();
            auto iter_corners();
        };

        auto iter_vertices()  const;
        auto iter_corners()   const;
        auto iter_halfedges() const;
        auto iter_facets()    const;
        auto iter_cells()     const;
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
        return cells.size() / nverts_per_cell();
    }

    inline int Volume::nfacets() const {
        return ncells() * nfacets_per_cell();
    }

    inline int Volume::nhalfedges() const {
        return ncells() * nhalfedges_per_cell();
    }

    inline constexpr int Volume::cell_from_facet(const int f) const {
        assert(f>=0 && f<nfacets());
        return f / nfacets_per_cell();
    }

    inline constexpr int Volume::cell_from_corner(const int c) const {
        assert(c>=0 && c<ncorners());
        return c / nverts_per_cell();
    }

    inline int Volume::vert(const int c, const int lv) const {
        assert(c>=0 && c<ncells() && lv>=0 && lv<nverts_per_cell());
        return cells[c*nverts_per_cell() + lv];
    }

    inline int &Volume::vert(const int c, const int lv) {
        assert(c>=0 && c<ncells() && lv>=0 && lv<nverts_per_cell());
        return cells[c*nverts_per_cell() + lv];
    }

    inline constexpr int Volume::nverts_per_cell() const {
        return reference_cells[cell_type].nverts();
    }

    inline constexpr int Volume::nfacets_per_cell() const {
        return reference_cells[cell_type].nfacets();
    }

    inline constexpr int Volume::nhalfedges_per_cell() const {
        return reference_cells[cell_type].ncorners();
    }

    inline constexpr int Volume::facet_size(const int f) const {
        assert(f>=0 && f<nfacets());
        return reference_cells[cell_type].facet_size(f % nfacets_per_cell());
    }

    inline int Volume::facet_vert(const int c, const int lf, const int lv) const {
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

    //////////////////////////////////////////////////////////////
    //                  _           _ _   _                     //
    //       _ __  _ __(_)_ __ ___ (_) |_(_)_   _____  ___      //
    //      | '_ \| '__| | '_ ` _ \| | __| \ \ / / _ \/ __|     //
    //      | |_) | |  | | | | | | | | |_| |\ V /  __/\__ \     //
    //      | .__/|_|  |_|_| |_| |_|_|\__|_| \_/ \___||___/     //
    //      |_|                                                 //
    //////////////////////////////////////////////////////////////

    inline Volume::Primitive::operator int() const {
        return id;
    }

    inline Volume::Primitive::operator int& () {
        return id;
    }

    inline bool Volume::Primitive::active() const {
        return id>=0;
    }

    inline Volume::Primitive& Volume::Primitive::operator=(const Volume::Primitive&& p) noexcept {
        return Primitive::operator=(p);
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

    inline vec3  Volume::Vertex::pos() const {
        return m.points[id];
    }

    inline vec3& Volume::Vertex::pos() {
        return const_cast<vec3 &>(m.points[id]); // TODO attention!
    }

    inline Volume::Vertex::operator const vec3&() const {
        return { m.points[id] };
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline Volume::Vertex Volume::Corner::vertex() const {
        return { m, m.vert(id / m.nverts_per_cell(), id % m.nverts_per_cell()) };
    }

    inline Volume::Cell Volume::Corner::cell() const {
         return { m, id / m.nverts_per_cell() };
    }

    inline int Volume::Corner::id_in_cell() const {
        return id % m.nverts_per_cell();
    }

    inline Volume::Halfedge Volume::Corner::halfedge() const {
        const int local_halfedge = reference_cells[m.cell_type].v2h[vertex() % m.nverts_per_cell()];
        return { m,  cell() * m.nhalfedges_per_cell() + local_halfedge };
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Volume::Halfedge::id_in_facet() const {
        return id - facet().halfedge(0);
    }

    inline int Volume::Halfedge::id_in_cell() const {
        return id % m.nhalfedges_per_cell();
    }

    inline Volume::Vertex Volume::Halfedge::from() const {
        return { m, m.vert(cell(), reference_cells[m.cell_type].from(id_in_cell())) };
    }

    inline Volume::Vertex Volume::Halfedge::to() const {
        return next().from();
    }

    inline Volume::Corner Volume::Halfedge::from_corner() const {
        return { m, cell() * m.nverts_per_cell() + reference_cells[m.cell_type].from(id_in_cell()) };
    }

    inline Volume::Corner Volume::Halfedge::to_corner() const {
        return next().from_corner();
    }

    inline Volume::Halfedge Volume::Halfedge::next() const {
        if (id_in_facet() < facet().size()-1)
            return { m, id + 1 };
        return { m, id - facet().size() + 1 };
    }

    inline Volume::Halfedge Volume::Halfedge::prev() const {
        if (id_in_facet())
            return { m, id - 1 };
        return { m, id + facet().size() - 1 };
    }

    inline Volume::Halfedge Volume::Halfedge::opposite_f() const {
        return { m, cell() * m.nhalfedges_per_cell() + reference_cells[m.cell_type].opposite(id_in_cell()) };
    }

    inline Volume::Halfedge Volume::Halfedge::opposite_c() const {
        assert(m.connected());
        Facet oppf = facet().opposite();
        if (!oppf.active()) return { m, -1 };
        for (int lv=0; lv<oppf.size(); lv++) {
            Halfedge res = oppf.halfedge(lv);
            if (res.from() == to() && res.to() == from()) return res;
        }
        um_assert(false);
        return { m, -1 };
    }

    inline Volume::Facet Volume::Halfedge::facet() const {
        return { m, cell() * m.nfacets_per_cell() + reference_cells[m.cell_type].facet(id_in_cell()) };
    }

    inline Volume::Cell Volume::Halfedge::cell() const {
        return { m, id / m.nhalfedges_per_cell() };
    }

    inline Volume::Halfedge::operator Segment3() const {
        return { from(), to() };
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
                if (!ref.active()) return iterator{ ref,ref };
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

    inline int Volume::Facet::size() const {
        return reference_cells[m.cell_type].facet_size(id % m.nfacets_per_cell());
    }

    inline int Volume::Facet::nverts() const {
        return m.facet(id).size();
    }

    inline int Volume::Facet::ncorners() const {
        return m.facet(id).size();
    }

    inline int Volume::Facet::nhalfedges() const {
        return m.facet(id).size();
    }

    inline Volume::Halfedge Volume::Facet::halfedge(int i) const{
        assert(i>=0 && i<size());
        return { m, cell() * m.nhalfedges_per_cell() + reference_cells[m.cell_type].corner(id_in_cell(), i) };
    }

    inline Volume::Vertex Volume::Facet::vertex(int i) const {
        assert(i>=0 && i<size());
        return { m, m.vert(cell(), reference_cells[m.cell_type].vert(id_in_cell(), i)) };
    }

    inline Volume::Corner Volume::Facet::corner(int i) const {
        assert(i>=0 && i<size());
        return { m, halfedge(i).from_corner() };
    }

    inline Volume::Facet Volume::Facet::opposite() const {
        assert(m.connected());
        return { m, m.conn->adjacent[id] };
    }

    inline bool Volume::Facet::on_boundary() const {
        return !opposite().active();
    }

    inline Volume::Cell Volume::Facet::cell() const {
        return { m, id / m.nfacets_per_cell() };
    }

    inline int Volume::Facet::id_in_cell() const {
        return id % m.nfacets_per_cell();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        return m.nhalfedges_per_cell();
    }

    inline Volume::Vertex Volume::Cell::vertex(int lv) const {
        assert(lv>=0 && lv<m.nverts_per_cell());
        return { m, m.vert(id, lv) };
    }

    inline Volume::Corner Volume::Cell::corner(int lc) const {
        assert(lc>=0 && lc<m.nverts_per_cell());
        return { m, id * m.nverts_per_cell() + lc };
    }

    inline Volume::Halfedge Volume::Cell::halfedge(int lh) const {
        assert(lh>=0 && lh<m.nhalfedges_per_cell());
        return { m, m.nhalfedges_per_cell()*id + lh };
    }

    inline Volume::Facet Volume::Cell::facet(int lf) const {
        assert(lf>=0 && lf<m.nfacets_per_cell());
        return { m, m.nfacets_per_cell()*id + lf };
    }

    inline Volume::Cell::operator Tetrahedron() const {
        um_assert(m.cell_type == Volume::TETRAHEDRON);
        return { vertex(0), vertex(1), vertex(2), vertex(3) };
    }

    inline Volume::Cell::operator Pyramid() const {
        um_assert(m.cell_type == Volume::PYRAMID);
        return { vertex(0), vertex(1), vertex(2), vertex(3), vertex(4) };
    }

    inline Volume::Cell::operator Hexahedron() const {
        um_assert(m.cell_type == Volume::HEXAHEDRON);
        return { vertex(0), vertex(1), vertex(2), vertex(3), vertex(4), vertex(5), vertex(6), vertex(7) };
    }

    inline Volume::Cell::operator Wedge() const {
        um_assert(m.cell_type == Volume::WEDGE);
        return { vertex(0), vertex(1), vertex(2), vertex(3), vertex(4), vertex(5) };
    }

    inline Volume::Facet::operator Triangle3() const {
        um_assert(size()==3);
        return { vertex(0), vertex(1), vertex(2) };
    }

    inline Volume::Facet::operator Quad3() const {
        um_assert(size()==4);
        return { vertex(0), vertex(1), vertex(2), vertex(3) };
    }

    inline Volume::Facet::operator Poly3() const {
        std::vector<vec3> pts(size());
        for (int i = 0; i < size(); i++)
            pts[i] = vertex(i);
        return { pts };
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    template <typename T> void Volume::delete_cells(const T &to_kill) {
        constexpr bool invocable = std::is_invocable_r_v<bool, T, int>;
        if constexpr (!invocable)
            assert(to_kill.size()==(size_t)ncells());

        std::vector<int>   cells_old2new(ncells(),   -1);
        std::vector<int>  facets_old2new(nfacets(),  -1);
        std::vector<int> corners_old2new(ncorners(), -1);

        int new_nb_cells   = 0;
        int new_nb_facets  = 0;
        int new_nb_corners = 0;

        for (auto c : iter_cells()) {
            if constexpr (invocable) {
                if (to_kill(c)) continue;
            } else {
                if (to_kill[c]) continue;
            }

            for (auto f : c.iter_facets())
                facets_old2new[f] = new_nb_facets++;
            for (auto corner : c.iter_corners()) {
                corners_old2new[corner] = new_nb_corners;
                cells[new_nb_corners] = corner.vertex();
                new_nb_corners++;
            }
            cells_old2new[c] = new_nb_cells++;
        }

        cells.resize(new_nb_corners);

        std::erase_if(attr_cells,   [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; }); // remove dead attributes
        std::erase_if(attr_facets,  [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; });
        std::erase_if(attr_corners, [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; });

        for (auto &wp : attr_cells)   if (auto spt = wp.lock())
            spt->compress(cells_old2new);
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->compress(facets_old2new);
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->compress(corners_old2new);
    }
}

#endif //__VOLUME_H__
