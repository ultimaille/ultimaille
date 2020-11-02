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

//      int  vert(const int c, const int lv) const;
//      int &vert(const int c, const int lv);
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

//      int  vert(const int c, const int lv) const;
//      int &vert(const int c, const int lv);
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

    inline int Volume::nverts() const {
        return points.size();
    }

    inline int Volume::vert(const int c, const int lv) const {
        return cells[corner(c, lv)];
    }

    inline int &Volume::vert(const int c, const int lv) {
        return cells[corner(c, lv)];
    }


}

#endif //__VOLUME_H__

