#ifndef __OBJ_H__
#define __OBJ_H__

#include <vector>
#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"
#include "ultimaille/attr_binding.h"

namespace UM {
    // TODO: export vn and vt attributes
    PolyLineAttributes read_wavefront_obj(const std::string &filename, PolyLine &m);
    SurfaceAttributes read_wavefront_obj(const std::string &filename, Triangles &m);
    SurfaceAttributes read_wavefront_obj(const std::string &filename, Quads  &m);
    SurfaceAttributes read_wavefront_obj(const std::string &filename, Polygons  &m);
    void write_wavefront_obj(const std::string &filename, const Surface &m, const SurfaceAttributes &attr = {});
    void write_wavefront_obj(const std::string &filename, const PolyLine &m, const PolyLineAttributes &attr = {});
}

#endif // __OBJ_H__

