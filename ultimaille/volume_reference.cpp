#include "volume_reference.h"

namespace UM {
    std::array<Polygons,4> reference_cells = {};
    std::array<SurfaceConnectivity,4> reference_conn = []() {
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

        return std::array<SurfaceConnectivity,4>{ tet, hex, wedge, pyramid };
    }();
}

