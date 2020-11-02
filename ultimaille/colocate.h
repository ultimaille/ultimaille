#ifndef __COLOCATE_H__
#define __COLOCATE_H__

#include <vector>
#include "geometry.h"

namespace UM {
    void colocate(const std::vector<vec3> &points, std::vector<int> &old2new, double tolerance);
}

#endif // __COLOCATE_H__

