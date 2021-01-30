#undef NDEBUG
#include <cassert>
#include <iostream>

#include <ultimaille/all.h>

using namespace UM;


int main() {
#if 1
    {
        Tetrahedra m;
        *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
        m.cells = {0,1,2,3};

        double vol = m.util.cell_volume(0);

        std::cerr << "Tet volume: " << vol << std::endl;
        assert(vol>0);

        for (int lf : range(4)) {
            vec3 n = m.util.facet_normal(0, lf);
            std::cerr << "Normal: " << n << std::endl;
//          assert(n*(m.util.bary_facet(0, lf) - m.util.bary_verts(0))>0);
        }
        write_geogram("tet.geogram", m, {{}, {}, {}, {}});
//      Tetrahedra z;
//      z = m;
//      std::cerr << z.ncells() << std::endl;
    }
#endif

    {
        Hexahedra m;
        *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1} };
        m.cells = {0,1,2,3,4,5,6,7};

        double vol = m.util.cell_volume(0);
        std::cerr << "Hex volume: " << vol << std::endl;
        //assert(vol>0);
        write_geogram("hex.geogram", m, {{}, {}, {}, {}});
    }
#if 0

    std::cerr << "Tet volume: " << vol << std::endl;
    for (int lf : range(4)) {
        vec3 n = facet_normal(m, 0, lf);
        assert(n*(facet_barycenter(m, 0, lf) - barycenter(m, 0))>0);
        std::cerr << n << std::endl;
    }
#endif


    return 0;
}

