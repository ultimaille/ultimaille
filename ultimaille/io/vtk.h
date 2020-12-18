#ifndef __VTK_H__
#define __VTK_H__

#include <vector>
#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"


// TODO : ATTRIBUTE NOT YET SUPPORTED

namespace UM {
    void write_vtk(const std::string filename, const PolyLine &pl);
    void write_vtk(const std::string filename, const Surface &m);
    // bool -> see medit.h
    void write_vtk(const std::string filename, const Volume  &m, bool hexes_GMSH_numerotation = true);


    void read_vtk(const std::string filename, PolyLine   &m);
    void read_vtk(const std::string filename, Triangles  &m);
    void read_vtk(const std::string filename, Quads      &m);
    void read_vtk(const std::string filename, Polygons   &m);
    void read_vtk(const std::string filename, Tetrahedra &m);
    void read_vtk(const std::string filename, Hexahedra  &m);
}


#endif // __MEDIT_H__

