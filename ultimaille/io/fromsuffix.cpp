#include "ultimaille/io/fromsuffix.h"
#include "ultimaille/io/geogram.h"
#include "ultimaille/io/medit.h"
#include "ultimaille/io/vtk.h"
#include "ultimaille/io/obj.h"



namespace UM {
    void write_fromsuffix(const std::string filename, const PolyLine& pl, PolyLineAttributes attr) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            write_medit(filename, pl);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            write_geogram(filename, pl, attr);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            write_vtk(filename, pl);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
    }
    void write_fromsuffix(const std::string filename, const Surface& m, SurfaceAttributes attr) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            write_medit(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            write_geogram(filename, m, attr);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            write_vtk(filename, m);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
    }
    void write_fromsuffix(const std::string filename, const Volume& m, VolumeAttributes attr) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            write_medit(filename, m);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            write_vtk(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            write_geogram(filename, m, attr);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
    }

    PolyLineAttributes read_fromsuffix(const std::string filename, PolyLine& m) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        PolyLineAttributes atts;
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            read_medit(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            atts = read_geogram(filename, m);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            read_vtk(filename, m);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
        return atts;
    }
    SurfaceAttributes  read_fromsuffix(const std::string filename, Triangles& m) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        SurfaceAttributes atts;
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            read_medit(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            atts = read_geogram(filename, m);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            read_vtk(filename, m);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
        return atts;
    }
    SurfaceAttributes  read_fromsuffix(const std::string filename, Quads& m) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        SurfaceAttributes atts;
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            read_medit(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            atts = read_geogram(filename, m);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            read_vtk(filename, m);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
        return atts;
    }
    SurfaceAttributes  read_fromsuffix(const std::string filename, Polygons& m) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        SurfaceAttributes atts;
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            read_medit(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            atts = read_geogram(filename, m);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            read_vtk(filename, m);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
        return atts;
    }
    VolumeAttributes   read_fromsuffix(const std::string filename, Tetrahedra& m) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        VolumeAttributes atts;
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            read_medit(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            atts = read_geogram(filename, m);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            read_vtk(filename, m);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
        return atts;
    }
    VolumeAttributes   read_fromsuffix(const std::string filename, Hexahedra& m) {
        if (filename.find(".") == std::string::npos) {
            std::cerr << "File name : \"" << filename << "\" has no suffix. " << std::endl;

        }
        VolumeAttributes atts;
        if (filename.size() > 5 && std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
            read_medit(filename, m);
        }
        else if (filename.size() > 8 && std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
            atts = read_geogram(filename, m);
        }
        else if (filename.size() > 4 && std::string(filename.end() - 4, filename.end()) == std::string(".vtk")) {
            read_vtk(filename, m);
        }
        else {
            std::cerr << "File format suffix not supported, use .geogram, .vtk, .obj or .mesh" << std::endl;
            exit(1);
        }
        return atts;
    }

}

