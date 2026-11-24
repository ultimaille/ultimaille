#include "primitive_geometry.h"

struct Volume::Connectivity {
    Volume &m;
    Volume::FacetAttribute<int> adjacent;

    Connectivity(Volume &m);
    void reset();
};

inline bool Volume::connected() const {
    return conn != nullptr;
}

struct Volume::Primitive {
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

struct Volume::Vertex : Primitive {
    using Primitive::Primitive;
    using Primitive::operator=;

    [[deprecated]] vec3  pos() const;
    [[deprecated]] vec3& pos();

    inline operator const vec3&() const;
};

struct Volume::Corner : Primitive {
    using Primitive::Primitive;
    using Primitive::operator=;

    int id_in_cell()    const;
    Vertex vertex()     const;
    Cell cell()         const;
    Halfedge halfedge() const;

    //          auto iter_halfedges() const; // TODO
};

struct Volume::Halfedge : Primitive {
    using Primitive::Primitive;
    using Primitive::operator=;

    int id_in_facet()     const;
    int id_in_cell()      const;
    Vertex from()         const;
    Vertex to()           const;
    Vertex vertex(int lv) const;
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

struct Volume::Facet : Primitive {
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

struct Volume::Cell : Primitive {
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

inline Volume::Vertex     Volume::vertex(int id) const { return   Vertex(*this, id); }
inline Volume::Corner     Volume::corner(int id) const { return   Corner(*this, id); }
inline Volume::Halfedge Volume::halfedge(int id) const { return Halfedge(*this, id); }
inline Volume::Facet       Volume::facet(int id) const { return    Facet(*this, id); }
inline Volume::Cell         Volume::cell(int id) const { return     Cell(*this, id); }

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
    return { m, m.vert(cell(), id_in_cell()) };
}

inline Volume::Cell Volume::Corner::cell() const {
    return { m, id / m.nverts_per_cell() };
}

inline int Volume::Corner::id_in_cell() const {
    return id % m.nverts_per_cell();
}

inline Volume::Halfedge Volume::Corner::halfedge() const {
    const int local_halfedge = reference_cells[m.cell_type].v2h[id_in_cell()];
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

inline Volume::Vertex Volume::Halfedge::vertex(int lv) const {
    assert(lv == 0 || lv == 1);
    if (lv) return to();
    return from();
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

