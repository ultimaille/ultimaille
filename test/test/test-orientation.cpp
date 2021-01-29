#undef NDEBUG
#include <cassert>
#include <iostream>

#include <ultimaille/all.h>

using namespace UM;


inline vec3 barycenter(const Volume &m, const int c) {
    vec3 ave = {0, 0, 0};
    const int nbv = m.nverts_per_cell();
    for (int lv=0; lv<nbv; lv++)
        ave += m.points[m.vert(c, lv)];
    return ave / static_cast<double>(nbv);
}

inline vec3 facet_barycenter(const Volume &m, const int c, const int lf) {
    vec3 ave = {0, 0, 0};
    const int nbv = m.facet_size(c, lf);
    for (int lv=0; lv<nbv; lv++)
        ave += m.points[m.facet_vert(c, lf, lv)];
    return ave / static_cast<double>(nbv);
}

// unit vector: weighted sum of normal of a triangle fan around the barycenter
inline vec3 facet_normal(const Volume &m, const int c, const int lf) {
    const int nbv = m.facet_size(c, lf);
    vec3 res = {0, 0, 0};
    vec3 bary = facet_barycenter(m, c, lf);
    for (int lv=0; lv<nbv; lv++)
        res += cross(
                m.points[m.facet_vert(c, lf,  lv       )]-bary,
                m.points[m.facet_vert(c, lf, (lv+1)%nbv)]-bary
                );
    return res.normalize();
}

int main() {
#if 1
     Tetrahedra m;
     *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
     m.cells = {0,1,2,3};

     double vol = m.cell_volume(0);
     std::cerr << "Tet volume: " << vol << std::endl;
     assert(vol>0);

     for (int lf : range(4)) {
         vec3 n = facet_normal(m, 0, lf);
         std::cerr << "Normal: " << n << std::endl;
         assert(n*(facet_barycenter(m, 0, lf) - barycenter(m, 0))>0);
     }
#endif

#if 0
     Hexahedra m;
     *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
     m.cells = {0,1,2,3};

     double vol = tet_volume(m, 0);
     assert(vol>0);

     std::cerr << "Tet volume: " << vol << std::endl;
     for (int lf : range(4)) {
         vec3 n = facet_normal(m, 0, lf);
         assert(n*(facet_barycenter(m, 0, lf) - barycenter(m, 0))>0);
         std::cerr << n << std::endl;
     }
#endif


    write_geogram("test.geogram", m, {{}, {}, {}, {}});
    return 0;
}

