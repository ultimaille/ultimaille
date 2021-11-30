#ifndef __VOLUME_REFERENCE_H__
#define __VOLUME_REFERENCE_H__

#include <array>
#include "surface.h"
#include "surface_connectivity.h"

namespace UM {
    extern std::array<Polygons,4> reference_cells;
    extern std::array<SurfaceConnectivity,4> reference_conn;
}

#endif //__VOLUME_REFERENCE_H__

