#ifndef __VOLUME_H__
#define __VOLUME_H__
#include <vector>
#include <memory>
#include "vec.h"
#include "pointset.h"

namespace UM {
    // signed volume for a tet with vertices (A,B,C,D) such that (AB, AC, AD) form a right hand basis
    inline double tet_volume(const vec3 &A, const vec3 &B, const vec3 &C, const vec3 &D) {
        return (cross(B-A, C-A)*(D-A))/6.;
    }

    struct GenericAttributeContainer;

    struct Volume { // polyhedral mesh interface
        PointSet points{};
        std::vector<int> cells{};

        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_cells{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_corners{};

        Volume() = default;

        int  create_cells(const int n);
        void delete_cells(const std::vector<bool> &to_kill);
        void delete_vertices(const std::vector<bool> &to_kill);

        void resize_attrs();
        void compress_attrs(const std::vector<bool> &cells_to_kill);

        int nverts()   const;
        int ncells()   const;
        int nfacets()  const;
        int ncorners() const;
        int cell_from_facet (const int f) const;
        int cell_from_corner(const int c) const;
        int  vert(const int c, const int lv) const;
        int &vert(const int c, const int lv);

        virtual int cell_type() const = 0; // TODO convert it to enum
        virtual int  nverts_per_cell() const = 0;
        virtual int nfacets_per_cell() const = 0;
//      virtual double cell_volume(const int c);

        virtual int facet_size(const int c, const int lf) const = 0;
        virtual int facet_vert(const int c, const int lf, const int lv) const = 0;
        virtual int  facet(const int c, const int lf) const = 0;
        virtual int corner(const int c, const int lc) const = 0;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Tetrahedra : Volume { // simplicial mesh implementation
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
        double cell_volume(const int c) const;
    };

    struct Hexahedra : Volume { // hex mesh implementation
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
    };

    struct Wedges : Volume { // prismatic mesh implementation
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct VolumeConnectivity { // half-edge-like connectivity interface
        VolumeConnectivity(const Volume &p_m);

        const Volume &m;
        std::vector<int> adjacent;
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

    inline int Volume::cell_from_facet (const int f) const {
        assert(f>=0 && f<nfacets());
        return f/nfacets_per_cell();
    }

    inline int Volume::cell_from_corner(const int c) const {
        assert(c>=0 && c<ncorners());
        return c/nverts_per_cell();
    }

    inline int Volume::vert(const int c, const int lv) const {
        return cells[corner(c, lv)];
    }

    inline int &Volume::vert(const int c, const int lv) {
        return cells[corner(c, lv)];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Tetrahedra::cell_type() const {
        return 0;
    }

    inline int Tetrahedra::nverts_per_cell() const {
        return 4;
    }

    inline int Tetrahedra::nfacets_per_cell() const {
        return 4;
    }

    inline int Tetrahedra::facet_size(const int c, const int lf) const {
        (void)c; (void)lf; // suppress unused parameter warnings
        return 3;
    }

    inline int Tetrahedra::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<4 && lv>=0 && lv<3);
        static constexpr int facet_vertex[4][3] = {{1,3,2}, {0,2,3}, {3,1,0}, {0,1,2}};
        return vert(c, facet_vertex[lf][lv]);
    }

    inline int Tetrahedra::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<4);
        return c*4 + lf;
    }

    inline int Tetrahedra::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<4);
        return c*4 + lc;
    }

    inline double Tetrahedra::cell_volume(const int c) const {
        return tet_volume(
                points[vert(c, 0)],
                points[vert(c, 1)],
                points[vert(c, 2)],
                points[vert(c, 3)]
                );
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Hexahedra::cell_type() const {
        return 1;
    }

    inline int Hexahedra::nverts_per_cell() const {
        return 8;
    }

    inline int Hexahedra::nfacets_per_cell() const {
        return 6;
    }

    inline int Hexahedra::facet_size(const int c, const int lf) const {
        (void)c; (void)lf; // suppress unused parameter warnings
        return 4;
    }

    inline int Hexahedra::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<6 && lv>=0 && lv<4);
        static constexpr int facet_vertex[6][4] = {{0,2,6,4}, {3,1,5,7}, {1,0,4,5}, {2,3,7,6}, {1,3,2,0}, {4,6,7,5}};
        return vert(c, facet_vertex[lf][lv]);
    }

    inline int Hexahedra::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<6);
        return c*6 + lf;
    }

    inline int Hexahedra::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<8);
        return c*8 + lc;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Wedges::cell_type() const {
        return 2;
    }

    inline int Wedges::nverts_per_cell() const {
        return 6;
    }

    inline int Wedges::nfacets_per_cell() const {
        return 5;
    }

    inline int Wedges::facet_size(const int c, const int lf) const {
        (void)c;
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        if (lf<2) return 3;
        return 4;
    }

    inline int Wedges::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell() && lv>=0 && lv<facet_size(c, lf));
        static constexpr int facet_vertex[5][4] = {{0,1,2,-1}, {3,5,4,-1}, {0,3,4,1}, {0,2,5,3}, {1,4,5,2} };
        return vert(c, facet_vertex[lf][lv]);
    }

    inline int Wedges::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        return c*6 + lf;
    }

    inline int Wedges::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<6);
        return c*6 + lc;
    }
}

#endif //__VOLUME_H__

