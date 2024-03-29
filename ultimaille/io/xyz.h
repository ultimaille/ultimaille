#ifndef __XYZ_H__
#define __XYZ_H__

#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/pointset.h"
#include "ultimaille/polyline.h"
#include "ultimaille/surface.h"
#include "ultimaille/volume.h"
#include "ultimaille/attr_binding.h"

namespace UM {
    void write_xyz(const std::string filename, const PointSet &ps);
    PointSetAttributes read_xyz(const std::string filename, PointSet &m);
}


#endif // __XYZ_H__

