#ifndef __MESH_IO_H__
#define __MESH_IO_H__

#include <vector>
#include <cstring>
#include "attributes.h"
#include "surface.h"
#include "polyline.h"

void read_wavefront_obj(const std::string filename, Polygons &m);
void write_wavefront_obj(const std::string filename, const Surface &m);

//void write_geogram_ascii(const std::string filename, const Surface &m, SurfaceAttributes attr = {{}, {}, {}});
void write_geogram(const std::string filename, const Surface &m, SurfaceAttributes attr = {{}, {}, {}});
void write_geogram(const std::string filename, const PolyLine &pl, PolyLineAttributes attr = {{}, {}});
SurfaceAttributes read_geogram(const std::string filename, Polygons &m);
PolyLineAttributes read_geogram(const std::string filename, PolyLine &pl);

#endif // __MESH_IO_H__

