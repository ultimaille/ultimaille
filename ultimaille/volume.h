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

    struct halfedge_around_edge_iter {
        const OppositeFacet  &of;
        const HalfEdgeHelper &heh;
        int start  = -1;
        int finish = -1;

        halfedge_around_edge_iter(const OppositeFacet &of, const int he);

        struct iterator {
            const OppositeFacet  &of;
            const HalfEdgeHelper &heh;
            int value;
            bool circ;

            void operator++() {
                value = of.opposite_c(heh.opposite_f(value));
                circ = false;
            }

            int operator*() const {
                return value;
            }

            bool operator!=(const iterator& rhs) const {
                return circ || value != rhs.value;
            }
        };

        iterator begin() const { return {of, heh, start, start==finish}; }
        iterator end()   const { return {of, heh, finish, false}; }
    };

    struct Volume {
        enum CELL_TYPE { TETRAHEDRON=0, HEXAHEDRON=1, WEDGE=2, PYRAMID=3 };
        constexpr virtual CELL_TYPE cell_type() const noexcept = 0;

        HalfEdgeHelper heh;
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

        Volume() : heh(*this), util(*this) {}
        Volume(const Volume& m) : heh(*this), util(*this) { // TODO re-think copying policy
            um_assert(!m.points.size() && !m.cells.size());
        }
        Volume& operator=(const Volume& m) {
            clear();
            um_assert(!m.points.size() && !m.cells.size());
            return *this;
        }

        struct Util {
            Util(const Volume &mesh) : m(mesh) {}
            virtual double cell_volume(const int c) const;
            virtual vec3 facet_normal(const int c, const int lf) const;
            vec3 bary_verts(const int c) const;
            vec3 bary_facet(const int c, const int lf) const;
            const Volume &m;
        } util;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Tetrahedra : Volume {
        constexpr CELL_TYPE cell_type() const noexcept override {
            return Volume::TETRAHEDRON;
        }
    };

    struct Hexahedra : Volume {
        constexpr CELL_TYPE cell_type() const noexcept override {
            return Volume::HEXAHEDRON;
        }
    };

    struct Wedges : Volume {
        constexpr CELL_TYPE cell_type() const noexcept override {
            return Volume::WEDGE;
        }
    };

    struct Pyramids : Volume {
        constexpr CELL_TYPE cell_type() const noexcept override {
            return Volume::PYRAMID;
        }
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
        return reference_cells[cell_type()].nverts();
    }

    inline constexpr int Volume::nfacets_per_cell() const {
        return reference_cells[cell_type()].nfacets();
    }

    inline constexpr int Volume::facet_size(const int f) const {
        assert(f>=0 && f<nfacets());
        return reference_cells[cell_type()].facet_size(f % nfacets_per_cell());
    }

    inline int Volume::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell() && lv>=0 && lv<facet_size(lf));
        return vert(c, reference_cells[cell_type()].vert(lf, lv));
    }

    inline constexpr int Volume::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        return c*nfacets_per_cell() + lf;
    }

    inline constexpr int Volume::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<nverts_per_cell());
        return c*nverts_per_cell() + lc;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline constexpr int HalfEdgeHelper::nhalfedges_per_cell() const {
        return reference_cells[m.cell_type()].ncorners();
    }

    // global halfedge id from local id
    inline constexpr int HalfEdgeHelper::halfedge(const int cell, const int cell_facet, const int facet_he) const {
//      assert(cell>=0 && cell<m.ncells());
        assert(cell_facet>=0 && cell_facet<m.nfacets_per_cell());
        assert(facet_he>=0 && facet_he<=m.facet_size(cell_facet));
        return cell*nhalfedges_per_cell() + reference_cells[m.cell_type()].corner(cell_facet, facet_he);
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
        return cell(he) * m.nverts_per_cell() + reference_cells[m.cell_type()].from(cell_halfedge(he));
    }

    // local facet id
    inline constexpr int HalfEdgeHelper::cell_facet(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return reference_cells[m.cell_type()].facet(cell_halfedge(he));
    }

    // local halfedge id
    inline constexpr int HalfEdgeHelper::cell_halfedge(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return he % nhalfedges_per_cell();
    }

    // local halfedge id
    inline constexpr int HalfEdgeHelper::facet_halfedge(const int he) const {
//      assert(he>=0 && he<nhalfedges());
        return cell_halfedge(he) - reference_cells[m.cell_type()].corner(cell_facet(he), 0);
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
        return nhalfedges_per_cell()*cell(he) + reference_cells[m.cell_type()].opposite(cell_halfedge(he));
    }


}

#endif //__VOLUME_H__

