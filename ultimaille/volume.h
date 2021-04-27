#ifndef __VOLUME_H__
#define __VOLUME_H__
#include <vector>
#include <memory>
#include "syntactic-sugar/assert.h"
#include "algebra/vec.h"
#include "pointset.h"

namespace UM {
    struct GenericAttributeContainer;

    struct Volume {
        enum CELL_TYPE { TETRAHEDRON=0, HEXAHEDRON=1, WEDGE=2, PYRAMID=3 };

        PointSet points{};
        std::vector<int> cells{};

        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_cells{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_corners{};

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

        virtual int cell_type() const = 0;
        virtual int  nverts_per_cell() const = 0;
        virtual int nfacets_per_cell() const = 0;

        virtual int facet_size(const int c, const int lf) const = 0;
        virtual int facet_vert(const int c, const int lf, const int lv) const = 0;
        virtual int  facet(const int c, const int lf) const = 0;
        virtual int corner(const int c, const int lc) const = 0;

        void clear() {
            points = {};
            cells  = {};
            attr_cells   = {};
            attr_facets  = {};
            attr_corners = {};
        }

        Volume() : util(*this) {}
        Volume(const Volume& m) : util(*this) {
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

    struct Tetrahedra : Volume {
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;

        Tetrahedra() : Volume(), util(*this) {}
        Tetrahedra(const Tetrahedra& m) : Volume(m), util(*this) {}
        Tetrahedra& operator=(const Tetrahedra& m) {
            Volume::operator=(m);
            return *this;
        }

        struct Util : Volume::Util {
            Util(const Tetrahedra &mesh) : Volume::Util(mesh) {}
            double cell_volume(const int c) const;
            vec3 facet_normal(const int c, const int lf) const;
        } util;
    };

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

    struct Hexahedra : Volume {
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
    };

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

    struct Wedges : Volume {
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
    };

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

    struct Pyramids : Volume {
        int cell_type() const;
        int  nverts_per_cell() const;
        int nfacets_per_cell() const;

        int facet_size(const int c, const int lf) const;
        int facet_vert(const int c, const int lf, const int lv) const;
        int  facet(const int c, const int lf) const;
        int corner(const int c, const int lc) const;
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
        return Volume::CELL_TYPE::TETRAHEDRON;
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
        static constexpr int facet_vertex[4][3] = {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}};
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
        return Volume::CELL_TYPE::HEXAHEDRON;
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
        assert(c>=0 && c<ncells()); assert(lf>=0 && lf<6); assert(lv>=0 && lv<4);
        static constexpr int facet_vertex[6][4] = {{0,4,6,2}, {1,3,7,5}, {0,1,5,4}, {2,6,7,3}, {0,2,3,1}, {4,5,7,6}};
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
        return Volume::CELL_TYPE::WEDGE;
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
        static constexpr int facet_vertex[5][4] = {{0,2,1,-1}, {3,4,5,-1}, {0,1,4,3}, {0,3,5,2}, {1,2,5,4}};
        return vert(c, facet_vertex[lf][lv]);
    }

    inline int Wedges::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        return c*5 + lf;
    }

    inline int Wedges::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<6);
        return c*6 + lc;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Pyramids::cell_type() const {
        return Volume::CELL_TYPE::PYRAMID;
    }

    inline int Pyramids::nverts_per_cell() const {
        return 5;
    }

    inline int Pyramids::nfacets_per_cell() const {
        return 5;
    }

    inline int Pyramids::facet_size(const int c, const int lf) const {
        (void)c;
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        if (!lf) return 4;
        return 3;
    }

    inline int Pyramids::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell() && lv>=0 && lv<facet_size(c, lf));
        static constexpr int facet_vertex[5][4] = {{0,3,2,1}, {0,1,4,-1}, {0,4,3,-1}, {2,3,4,-1}, {1,2,4,-1}};
        return vert(c, facet_vertex[lf][lv]);
    }

    inline int Pyramids::facet(const int c, const int lf) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell());
        return c*5 + lf;
    }

    inline int Pyramids::corner(const int c, const int lc) const {
        assert(c>=0 && c<ncells() && lc>=0 && lc<6);
        return c*5 + lc;
    }
}

#endif //__VOLUME_H__

