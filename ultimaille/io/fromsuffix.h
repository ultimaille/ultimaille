#ifndef __FROMSUFFIX_H__
#define __FROMSUFFIX_H__

#include <vector>
#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"



// TODO : ATTRIBUTE NOT YET SUPPORTED

namespace UM {
    void write_fromsuffix(const std::string filename, const PolyLine& pl, PolyLineAttributes attr = { {}, {} });
    void write_fromsuffix(const std::string filename, const Surface& m, SurfaceAttributes attr = { {}, {}, {} });
    void write_fromsuffix(const std::string filename, const Volume& m, VolumeAttributes attr = { {}, {}, {}, {} });

    PolyLineAttributes read_fromsuffix(const std::string filename, PolyLine& m);
    SurfaceAttributes  read_fromsuffix(const std::string filename, Triangles& m);
    SurfaceAttributes  read_fromsuffix(const std::string filename, Quads& m);
    SurfaceAttributes  read_fromsuffix(const std::string filename, Polygons& m);
    VolumeAttributes   read_fromsuffix(const std::string filename, Tetrahedra& m);
    VolumeAttributes   read_fromsuffix(const std::string filename, Hexahedra& m);

}


#endif // __MEDIT_H__

