#ifndef __VOLUME_H__
#define __VOLUME_H__
#include <vector>
#include <memory>
#include "geometry.h"
#include "pointset.h"

namespace UM {
    struct GenericAttributeContainer;

    struct Volume { // polyhedral mesh interface
        PointSet points{};
        std::vector<int> cells{};

        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_cells{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_corners{};

        Volume() = default;

        ////	void delete_vertices(const std::vector<bool> &to_kill);
        ////	virtual void delete_facets(const std::vector<bool> &to_kill);
        void resize_attrs();
        void compress_attrs(const std::vector<bool> &cells_to_kill);

        int nverts() const;
        int  vert(const int c, const int lv) const;
        int &vert(const int c, const int lv);

        virtual int ncells()   const = 0;
        virtual int nfacets()  const = 0;
        virtual int ncorners() const = 0;

        virtual int  nverts_per_cell() const = 0;
        virtual int nfacets_per_cell() const = 0;
        virtual int cell_type() const = 0;

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

        int ncorners() const;
        int ncells()  const;
        int nfacets() const;
        int create_tets(const int n);

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
    };

    struct Hexahedra : Volume { // hex mesh implementation
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int ncorners() const;
        int ncells()  const;
        int nfacets() const;
        int create_hexa(const int n);

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
    };

    struct VolumeConnectivity { // half-edge-like connectivity interface
        VolumeConnectivity(const Volume &p_m);

        const Volume &m;
        std::vector<int> adjacent;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Volume::nverts() const {
        return points.size();
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

    inline int Tetrahedra::ncorners() const {
        return ncells()*4;
    }

    inline int Tetrahedra::ncells() const {
        return cells.size()/4;
    }

    inline int Tetrahedra::nfacets() const {
        return ncells()*4;
    }

    inline int Tetrahedra::facet_size(const int c, const int lf) const {
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

    inline int Hexahedra::ncorners() const {
        return ncells()*8;
    }

    inline int Hexahedra::ncells() const {
        return cells.size()/8;
    }

    inline int Hexahedra::nfacets() const {
        return ncells()*6;
    }

    inline int Hexahedra::facet_size(const int c, const int lf) const {
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
}

#endif //__VOLUME_H__

