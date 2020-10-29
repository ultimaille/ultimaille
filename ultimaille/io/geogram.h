#ifndef __MESH_IO_H__
#define __MESH_IO_H__

#include <vector>
#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"

namespace UM {
    void read_wavefront_obj(const std::string filename, Polygons &m);
    void write_wavefront_obj(const std::string filename, const Surface &m);

    //void write_geogram_ascii(const std::string filename, const Surface &m, SurfaceAttributes attr = {{}, {}, {}});

    void write_geogram(const std::string filename, const PolyLine &pl, PolyLineAttributes attr = {{}, {}});
    void write_geogram(const std::string filename, const Surface &m, SurfaceAttributes attr = {{}, {}, {}});
    void write_geogram(const std::string filename, const Volume  &m,  VolumeAttributes attr = {{}, {}, {}, {}});

    PolyLineAttributes read_geogram(const std::string filename, PolyLine &m);
    SurfaceAttributes  read_geogram(const std::string filename, Polygons &m);
    VolumeAttributes   read_geogram(const std::string filename, Tetrahedra &m);
    VolumeAttributes   read_geogram(const std::string filename, Hexahedra  &m);
}

#endif // __MESH_IO_H__

