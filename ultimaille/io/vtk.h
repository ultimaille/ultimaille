#ifndef __VTK_H__
#define __VTK_H__

#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"

// TODO: colors and scalar fields are yet to be implemented

namespace UM {
/*
    void write_vtk(const std::string filename, const PointSet &ps, const PointSetAttributes attr = {});
    void write_vtk(const std::string filename, const PolyLine &pl, const PolyLineAttributes attr = {});
    void write_vtk(const std::string filename, const Triangles &m, const SurfaceAttributes attr = {});
    void write_vtk(const std::string filename, const Quads &m, const SurfaceAttributes attr = {});
    void write_vtk(const std::string filename, const Polygons &m, const SurfaceAttributes attr = {});
    void write_vtk(const std::string filename, const Tetrahedra &m, const VolumeAttributes attr = {});
    void write_vtk(const std::string filename, const Hexahedra  &m, const VolumeAttributes attr = {});
    void write_vtk(const std::string filename, const Wedges     &m, const VolumeAttributes attr = {});
    void write_vtk(const std::string filename, const Pyramids   &m, const VolumeAttributes attr = {});
*/
    void write_vtk(const std::string filename, const PointSet &ps, const PointSetAttributes attr = {});
    void write_vtk(const std::string filename, const PolyLine &pl, const PolyLineAttributes attr = {});
    void write_vtk(const std::string filename, const Surface &m,   const SurfaceAttributes  attr = {});
    void write_vtk(const std::string filename, const Volume  &m,   const VolumeAttributes   attr = {});


    PointSetAttributes read_vtk(const std::string filename, PointSet   &m);
    PolyLineAttributes read_vtk(const std::string filename, PolyLine   &m);
    SurfaceAttributes  read_vtk(const std::string filename, Triangles  &m);
    SurfaceAttributes  read_vtk(const std::string filename, Quads      &m);
    SurfaceAttributes  read_vtk(const std::string filename, Polygons   &m);
    VolumeAttributes   read_vtk(const std::string filename, Tetrahedra &m);
    VolumeAttributes   read_vtk(const std::string filename, Hexahedra  &m);
    VolumeAttributes   read_vtk(const std::string filename, Wedges     &m);
    VolumeAttributes   read_vtk(const std::string filename, Pyramids   &m);
}


#endif // __VTK_H__

