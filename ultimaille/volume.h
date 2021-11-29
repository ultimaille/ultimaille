#ifndef __VOLUME_H__
#define __VOLUME_H__
#include <vector>
#include <array>
#include <memory>
#include "syntactic-sugar/assert.h"
#include "algebra/vec.h"
#include "pointset.h"
#include "surface.h"
#include "surface_connectivity.h"

namespace UM {
    static std::array<Polygons,4> reference_cells = {};
    static auto reference_fec = []() -> std::array<SurfaceConnectivity,4> {
        /**
         *  LOCAL NUMBERING CONVENTION
         *
         *     Z             3            The vectors (01), (02), (03) form a right-hand basis.
         *     ^            /.\           The facets are numbered w.r.t the opposite vertex.
         *     |           / . \          For the reference tetrahedron v0=(0,0,0), v1=(1,0,0), v2=(0,1,0), v3=(0,0,1),
         *     |          /  .  \         the facets are numbered as follows:
         *     o         /   .   \        f0:x+y+z=1 f1:x=0 f2:y=0, f3:z=0
         *    / \       /    .    \
         *  /     \    /   . 0 .   \        The vertices inside each facet are numbered in a way that the normal vector
         * X       Y  /  .       .  \       points outside (CCW ordering when viewed from outside).
         *           / .           . \      The smallest local index is the first facet vertex.
         *          /_________________\
         *         1                   2
         */

        Polygons& tet = reference_cells[0];
        *tet.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
        tet.facets = {1,2,3, 0,3,2, 0,1,3, 0,2,1};
        tet.offset = {0,3,6,9,12};

        /**
         * LOCAL NUMBERING CONVENTION
         * The numbering of the vertices is derived from the binary         The facets are numbered in the following order:
         * code corresponding to the coordinates of the unit cube:          f0:x=0  f1:x=1
         * v0:(0,0,0), v1:(1,0,0), v2:(0,1,0), v3:(1,1,0)                   f2:y=0  f3:y=1
         * v4:(0,0,1), v5:(1,0,1), v6:(0,1,1), v7:(1,1,1)                   f4:z=0  f5:z=1
         *
         *           6-----------7                                                    +-----------+
         *          /|          /|                                                   /|          /|
         *         / |         / |                                                  / |   5     / |
         *        /  |        /  |                                                 /  |     3  /  |
         *       4-----------5   |                                                +-----------+   |
         *       |   |       |   |                                                | 0 |       | 1 |
         *       |   2-------|---3                                                |   +-------|---+
         *       |  /        |  /                                                 |  /  2     |  /
         * Z     | /         | /                                            Z     | /     4   | /
         * ^  Y  |/          |/                                             ^  Y  |/          |/
         * | /   0-----------1                                              | /   +-----------+
         * |/                                                               |/
         * o----> X                                                         o----> X
         *
         * The vertices inside each facet are numbered in a way that the normal vector points outside
         * (CCW ordering when viewed from outside).
         * The smallest local index is the first facet vertex.
         */

        Polygons& hex = reference_cells[1];
        *hex.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
        hex.facets = {0,4,6,2, 1,3,7,5, 0,1,5,4, 2,6,7,3, 0,2,3,1, 4,5,7,6};
        hex.offset = {0,4,8,12,16,20,24};

        /**
         * LOCAL NUMBERING CONVENTION
         *
         *            5                                         ^
         *           /. \                                      /. \
         *          / .   \                                   / .   \
         *         /  .     \                                /  .1    \
         *        /   .       \                             /   .       \
         *       3-------------4                           +-------------+
         *       |    .        |                           |  3 .   4    |
         *       |    .        |                           |    .        |
         *       |    2        |                           |    +        |
         *       |   .  .      |                           |   .  2      |
         *       |  .     .    |                           |  .     .    |
         * Z     | .        .  |                     Z     | .   0    .  |
         * ^  Y  |.           .|                     ^  Y  |.           .|
         * | /   0-------------1                     | /   +-------------+
         * |/                                        |/
         * o----> X                                  o----> X
         *
         * The vertices inside each facet are numbered in a way that the normal vector points outside
         * (CCW ordering when viewed from outside).
         * The smallest local index is the first facet vertex.
         */

        Polygons& wedge = reference_cells[2];
        *wedge.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {0,1,1}};
        wedge.facets = {0,2,1, 3,4,5, 0,1,4,3, 0,3,5,2, 1,2,5,4};
        wedge.offset = {0,3,6,10,14,18};

        /**
         * LOCAL NUMBERING CONVENTION
         *
         *                   4                                                      .
         *                 .::'.                                                  .::'.
         *                : : : '.                                               : : : '.
         *              .' :  :  '.                                            .' :  :  '.
         *             .'  :  :    '.                                         .'  :  :    '.
         *            :   :   :     '.                                       :   :   :     '.
         *          .:    :    :      :                                    .:    :    :      :
         *         .'    :     :       '.                                 .'    :3    :       '.
         *        .'     :     :         :                               .'     :     :    4    :
         *       :      :      :          '.                            :      :      :          '.
         *     .'       :       :          '.                         .'   2   :      :1          '.
         *    .'       :     .. 2            '.                      .'       :     ..'.            '.
         *   :  ......':'''''     '''...      '.                    :  ......':'''''     '''...      '.
         *  3.'''     :                 '''...  '                  ..'''     :                 '''...  '
         *   '.       :                       ''' 1                 '.       :      0                ''':
         *    '.     :                     ...'''           Z        '.     :                     ...'''
         *     '.    :               ...'''          Y     :          '.    :               ...'''
         *      '.  :         ...''''                 '.  :            '.  :         ...''''
         *       '. :   ...'''                         '. :   ..>X      '. :   ...'''
         *         `0.''                                 'o.''            `..''
         *
         *
         * The vertices inside each facet are numbered in a way that the normal vector points outside
         * (CCW ordering when viewed from outside).
         * The smallest local index is the first facet vertex.
         */

        Polygons& pyramid = reference_cells[3];
        *pyramid.points.data = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {.5,.5,.5}};
        pyramid.facets = {0,3,2,1, 0,1,4, 0,4,3, 2,3,4, 1,2,4};
        pyramid.offset = {0,4,7,10,13,16};

        return { tet, hex, wedge, pyramid };
    }();

    struct GenericAttributeContainer;

    struct Volume {
        enum CELL_TYPE { TETRAHEDRON=0, HEXAHEDRON=1, WEDGE=2, PYRAMID=3 };

        const CELL_TYPE cell_type;
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
        int cell_from_facet (const int f) const;
        int cell_from_corner(const int c) const;
        int  vert(const int c, const int lv) const;
        int &vert(const int c, const int lv);

        int  nverts_per_cell() const;
        int nfacets_per_cell() const;
        int facet_size(const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;

        void clear() {
            points = {};
            cells  = {};
            attr_cells   = {};
            attr_facets  = {};
            attr_corners = {};
        }

        Volume(CELL_TYPE cell_type) : cell_type(cell_type), util(*this) {}
        Volume(const Volume& m) : cell_type(m.cell_type), util(*this) { // TODO re-think copying policy
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

    inline int Volume::cell_from_facet(const int f) const {
        assert(f>=0 && f<nfacets());
        return f/nfacets_per_cell();
    }

    inline int Volume::cell_from_corner(const int c) const {
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

    inline int Volume::nverts_per_cell() const {
        return reference_cells[cell_type].nverts();
    }

    inline int Volume::nfacets_per_cell() const {
        return reference_cells[cell_type].nfacets();
    }

    inline int Volume::facet_size(const int lf) const {
        assert(lf>=0 && lf<nfacets_per_cell());
        return reference_cells[cell_type].facet_size(lf);
    }

    inline int Volume::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell() && lv>=0 && lv<facet_size(lf));
        return vert(c, reference_cells[cell_type].vert(lf, lv));
    }

    inline int Volume::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        return c*nfacets_per_cell() + lf;
    }

    inline int Volume::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<nverts_per_cell());
        return c*nverts_per_cell() + lc;
    }
}

#endif //__VOLUME_H__

