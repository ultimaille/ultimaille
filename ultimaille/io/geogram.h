#ifndef __GEOGRAM_H__
#define __GEOGRAM_H__

#include <vector>
#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/surface_connectivity.h"
#include "ultimaille/polyline.h"

namespace UM {
    void write_geogram(const std::string filename, const PointSet &ps, const PointSetAttributes attr = {{}});
    void write_geogram(const std::string filename, const PolyLine &pl, const PolyLineAttributes attr = {{}, {}});
    void write_geogram(const std::string filename, const Surface &m,   const SurfaceAttributes  attr = {{}, {}, {}});
    void write_geogram(const std::string filename, const Volume  &m,   const VolumeAttributes   attr = {{}, {}, {}, {}});

    PointSetAttributes read_geogram(const std::string filename, PointSet   &m);
    PolyLineAttributes read_geogram(const std::string filename, PolyLine   &m);
    SurfaceAttributes  read_geogram(const std::string filename, Triangles  &m);
    SurfaceAttributes  read_geogram(const std::string filename, Quads      &m);
    SurfaceAttributes  read_geogram(const std::string filename, Polygons   &m);
    VolumeAttributes   read_geogram(const std::string filename, Tetrahedra &m);
    VolumeAttributes   read_geogram(const std::string filename, Hexahedra  &m);
    VolumeAttributes   read_geogram(const std::string filename, Wedges     &m);
    VolumeAttributes   read_geogram(const std::string filename, Pyramids   &m);
}

#endif // __GEOGRAM_H__

