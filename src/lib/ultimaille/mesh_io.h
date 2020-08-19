#ifndef __MESH_IO_H__
#define __MESH_IO_H__

#include <vector>
#include <cstring>
#include "surface.h"

void read_wavefront_obj(const std::string filename, PolyMesh &m);
void write_wavefront_obj(const std::string filename, const Surface &m);
void write_geogram_ascii(const std::string filename, const Surface &m,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > fattr,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > cattr);

#endif // __MESH_IO_H__

