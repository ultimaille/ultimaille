#ifndef __VOLUME_REFERENCE_H__
#define __VOLUME_REFERENCE_H__
#include <array>

namespace UM {
    struct ReferenceCell {     // There are 4 types of volume meshes: Tetrahedra, Hexahedra, Wedges, Prisms.
        const int nv, nf, nc;  // Each one has a reference cell (a polygonal surface) that encodes the numbering convention.
        const vec3 points[8];  // See comments below about the conevntion.
        const int facets[24];  // This struct is a minimal polygonal surface interface needed to encode connectivity between cells.
        const int offset[7];
        const int c2f[24];
        const int opp[24];

        constexpr int nverts()   const { return nv; }
        constexpr int nfacets()  const { return nf; }
        constexpr int ncorners() const { return nc; }
        constexpr int facet_size(const int f) const { return offset[f+1]-offset[f]; }

        constexpr int corner(const int f, const int lc) const { return offset[f]+lc; }
        constexpr int vert(const int f, const int lv)   const { return facets[offset[f]+lv]; }

        constexpr int facet(const int he)    const { return c2f[he]; }
        constexpr int opposite(const int he) const { return opp[he]; }
        constexpr int from(const int he)     const { return vert(c2f[he], he - offset[c2f[he]]); }
    };

    constexpr std::array<ReferenceCell,4> reference_cells = {{
        {4,4,12,{{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{1,2,3,0,3,2,0,1,3,0,2,1},{0,3,6,9,12},{0,0,0,1,1,1,2,2,2,3,3,3},{10,4,7,8,1,9,11,2,3,5,0,6}},
        {8,6,24,{{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}},{0,4,6,2,1,3,7,5,0,1,5,4,2,6,7,3,0,2,3,1,4,5,7,6},{0,4,8,12,16,20,24},{0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5},{11,23,12,16,18,14,21,9,19,7,20,0,2,22,5,17,3,15,4,8,10,6,13,1}},
        {6,5,18,{{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}},{0,2,1,3,4,5,0,1,4,3,0,3,5,2,1,2,5,4},{0,3,6,10,14,18},{0,0,0,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4},{13,14,6,8,16,11,2,17,3,10,9,5,15,0,1,12,4,7}},
        {5,5,16,{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0.5,0.5,0.5}},{0,3,2,1,0,1,4,0,4,3,2,3,4,1,2,4},{0,4,7,10,13,16},{0,0,0,0,1,1,1,2,2,2,3,3,3,4,4,4},{9,10,13,4,3,15,7,6,11,0,1,8,14,2,12,5}}
    }};
}

#endif //__VOLUME_REFERENCE_H__

