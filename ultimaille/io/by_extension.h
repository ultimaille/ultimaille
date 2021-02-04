#ifndef __BY_EXTENSION_H__
#define __BY_EXTENSION_H__

#include <filesystem>
#include <type_traits>

#include "ultimaille/attributes.h"
#include "ultimaille/io/geogram.h"
#include "ultimaille/io/medit.h"
#include "ultimaille/io/vtk.h"
#include "ultimaille/io/obj.h"

namespace UM {
    inline PolyLineAttributes empty_attr(const PolyLine &m) {
        (void)m;
        return {};
    }

    inline SurfaceAttributes empty_attr(const Surface &m) {
        (void)m;
        return {};
    }

    inline VolumeAttributes empty_attr(const Volume &m) {
        (void)m;
        return {};
    }

    template <class M> void write_by_extension(const std::string path, const M &m, const decltype(empty_attr(m)) a = {}) {
        std::string ext = std::filesystem::path(path).extension().string();
        if (ext == ".geogram")
            write_geogram(path, m, a);
        if (ext == ".mesh")
            write_medit(path, m);
        if (ext == ".vtk")
            write_vtk(path, m);
        if constexpr (std::is_same_v<decltype(empty_attr(m)), SurfaceAttributes>) {
            if (ext == ".obj")
                write_wavefront_obj(path, m, a);
        }
    }

    template <class M> auto read_by_extension(const std::string path, M& m) -> decltype(empty_attr(m)) {
        std::string ext = std::filesystem::path(path).extension().string();
        if (ext == ".geogram")
            return read_geogram(path, m);
        if (ext == ".mesh")
            return read_medit(path, m);
        if (ext == ".vtk")
            return read_vtk(path, m);
        if constexpr (std::is_same_v<decltype(empty_attr(m)), SurfaceAttributes>) {
            if (ext == ".obj")
                return read_wavefront_obj(path, m);
        }
        return {};
    }
}

#endif // __BY_EXTENSION_H__

