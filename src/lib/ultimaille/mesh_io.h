#ifndef __MESH_IO_H__
#define __MESH_IO_H__

#include <vector>
#include <cstring>
#include "surface.h"

void read_wavefront_obj(const std::string filename, PolyMesh &m);
void write_wavefront_obj(const std::string filename, const Surface &m);

typedef std::tuple<std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > >,
                   std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > >,
                   std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > > SurfaceAttributes;

void write_geogram_ascii(const std::string filename, const Surface &m, SurfaceAttributes attr = {{}, {}, {}});
void write_geogram(const std::string filename, const Surface &m, SurfaceAttributes attr = {{}, {}, {}});
SurfaceAttributes read_geogram(const std::string filename, PolyMesh &m);

#endif // __MESH_IO_H__

