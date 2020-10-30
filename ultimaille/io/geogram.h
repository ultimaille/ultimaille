#ifndef __GEOGRAM_H__
#define __GEOGRAM_H__

#include <vector>
#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"

namespace UM {
    void write_geogram(const std::string filename, const PolyLine &pl, PolyLineAttributes attr = {{}, {}});
    void write_geogram(const std::string filename, const Surface &m, SurfaceAttributes attr = {{}, {}, {}});
    void write_geogram(const std::string filename, const Volume  &m,  VolumeAttributes attr = {{}, {}, {}, {}});

    PolyLineAttributes read_geogram(const std::string filename, PolyLine   &m);
    SurfaceAttributes  read_geogram(const std::string filename, Triangles  &m);
    SurfaceAttributes  read_geogram(const std::string filename, Quads      &m);
    SurfaceAttributes  read_geogram(const std::string filename, Polygons   &m);
    VolumeAttributes   read_geogram(const std::string filename, Tetrahedra &m);
    VolumeAttributes   read_geogram(const std::string filename, Hexahedra  &m);
}

#endif // __GEOGRAM_H__

