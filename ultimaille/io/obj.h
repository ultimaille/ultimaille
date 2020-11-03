#ifndef __OBJ_H__
#define __OBJ_H__

#include <vector>
#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"

namespace UM {
    // TODO: export vn and vt attributes
    void read_wavefront_obj(const std::string filename, Triangles &m);
    void read_wavefront_obj(const std::string filename, Polygons  &m);
    void write_wavefront_obj(const std::string filename, const Surface &m);
}

#endif // __OBJ_H__

