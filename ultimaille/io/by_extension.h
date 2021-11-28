#ifndef __BY_EXTENSION_H__
#define __BY_EXTENSION_H__

#include <filesystem>
#include <type_traits>

#include "ultimaille/attributes.h"
#include "ultimaille/io/geogram.h"
#include "ultimaille/io/medit.h"
#include "ultimaille/io/vtk.h"
#include "ultimaille/io/xyz.h"
#include "ultimaille/io/obj.h"

namespace UM {
    inline PointSetAttributes empty_attr([[maybe_unused]] const PointSet &m) { return {}; }
    inline PolyLineAttributes empty_attr([[maybe_unused]] const PolyLine &m) { return {}; }
    inline SurfaceAttributes  empty_attr([[maybe_unused]] const Surface  &m) { return {}; }
    inline VolumeAttributes   empty_attr([[maybe_unused]] const Volume   &m) { return {}; }

    template <class M> void write_by_extension(const std::string path, const M &m, const decltype(empty_attr(m)) a = {}) {
        std::string ext = std::filesystem::path(path).extension().string();
        try {
            if (ext == ".geogram")
                write_geogram(path, m, a);
            if (ext == ".mesh")
                write_medit(path, m);
            if (ext == ".vtk")
                write_vtk(path, m, a);
            if constexpr (std::is_same_v<decltype(empty_attr(m)), PointSetAttributes>) {
                if (ext == ".xyz")
                    write_xyz(path, m);
            }
            if constexpr (std::is_same_v<decltype(empty_attr(m)), SurfaceAttributes>) {
                if (ext == ".obj")
                    write_wavefront_obj(path, m, a);
            }
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    template <class M> auto read_by_extension(const std::string path, M& m) -> decltype(empty_attr(m)) {
        std::string ext = std::filesystem::path(path).extension().string();
        try {
            if (ext == ".geogram")
                return read_geogram(path, m);
            if (ext == ".mesh")
                return read_medit(path, m);
            if (ext == ".vtk")
                return read_vtk(path, m);
            if constexpr (std::is_same_v<decltype(empty_attr(m)), PointSetAttributes>) {
                if (ext == ".xyz")
                    return read_xyz(path, m);
            }
            if constexpr (std::is_same_v<decltype(empty_attr(m)), SurfaceAttributes>) {
                if (ext == ".obj")
                    return read_wavefront_obj(path, m);
            }
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
        }
        return {};
    }
}

#endif // __BY_EXTENSION_H__

