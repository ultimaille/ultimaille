#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "pointset.h"
#include "syntactic-sugar/assert.h"

namespace UM {
    struct GenericAttributeContainer;

    struct Surface { // polygonal mesh interface
        PointSet points{};
        std::vector<int> facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_facets{};
        std::vector<std::weak_ptr<GenericAttributeContainer> > attr_corners{};

        void delete_vertices(const std::vector<bool> &to_kill);
        virtual void delete_facets(const std::vector<bool> &to_kill);
        void delete_isolated_vertices();
        void resize_attrs();
        void compress_attrs(const std::vector<bool> &facets_to_kill);

        int nverts() const;
        int ncorners() const;

        virtual int nfacets() const = 0;
        virtual int facet_size(const int fi) const = 0;
        virtual int corner(const int fi, const int ci) const = 0;
        virtual int   vert(const int fi, const int lv) const = 0;
        virtual int  &vert(const int fi, const int lv)       = 0;

        virtual void clear() {
            points = {};
            attr_facets  = {};
            attr_corners = {};
        }

        Surface() : util(*this) {}
        Surface(const Surface& m) : util(*this) {
            um_assert(!m.points.size() && !m.facets.size());
        }
        Surface& operator=(const Surface& m) {
            clear();
            um_assert(!m.points.size() && !m.facets.size());
            return *this;
        }

        struct Util {
            Util(const Surface &mesh) : m(mesh) {}
            vec3 bary_verts(const int f) const;
            const Surface &m;
        } util;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Triangles : Surface { // simplicial mesh implementation
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        Triangles() : Surface(), util(*this) {}
        Triangles(const Triangles& m) : Surface(m), util(*this) {}
        Triangles& operator=(const Triangles& m) {
            Surface::operator=(m);
            return *this;
        }

        struct Util : Surface::Util {
            Util(const Triangles &mesh) : Surface::Util(mesh) {}
            double unsigned_area(const int f) const;
            void project(const int t, vec2& z0, vec2& z1, vec2& z2) const;
            vec3 normal(const int f) const;
        } util;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Quads : Surface { // quad mesh implementation
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        Quads() : Surface(), util(*this) {}
        Quads(const Quads& m) : Surface(m), util(*this) {}
        Quads& operator=(const Quads& m) {
            Surface::operator=(m);
            return *this;
        }

        struct Util : Surface::Util {
            Util(const Quads &mesh) : Surface::Util(mesh) {}
            double unsigned_area(const int f) const;
            vec3 normal(const int f) const;
        } util;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Polygons : Surface { // polygonal mesh implementation
        std::vector<int> offset = {};
        Polygons() : Surface(), offset(1, 0), util(*this)  {}
        Polygons(const Polygons& m) : Surface(m), offset(1,0), util(*this) {}
        Polygons& operator=(const Polygons& m) {
            clear();
            um_assert(!m.points.size() && !m.facets.size() && 1==m.offset.size());
            return *this;
        }

        int create_facets(const int n, const int size);
        void delete_facets(const std::vector<bool> &to_kill);

        virtual void clear() {
            Surface::clear();
            offset = {0};
        }

        int nfacets()  const;
        int facet_size(const int fi) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        struct Util : Surface::Util {
            Util(const Polygons &mesh) : Surface::Util(mesh) {}
        } util;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // these implementations are here and not in the .cpp because all inline functions must be available in all translation units //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Surface::nverts() const {
        return points.size();
    }

    inline int Surface::ncorners() const {
        return facets.size();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Triangles::nfacets() const {
        assert(0==facets.size()%3);
        return facets.size()/3;
    }

    inline int Triangles::facet_size(const int) const {
        return 3;
    }

    inline int Triangles::corner(const int fi, const int ci) const {
        assert(ci>=0 && ci<3 && fi>=0 && fi<nfacets());
        return fi*3 + ci;
    }

    inline int Triangles::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
        return facets[fi*3 + lv];
    }

    inline int &Triangles::vert(const int fi, const int lv) {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<3);
        return facets[fi*3 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Quads::nfacets() const {
        assert(0==facets.size()%4);
        return facets.size()/4;
    }

    inline int Quads::facet_size(const int) const {
        return 4;
    }

    inline int Quads::corner(const int fi, const int ci) const {
        assert(ci>=0 && ci<4 && fi>=0 && fi<nfacets());
        return fi*4 + ci;
    }

    inline int Quads::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<4);
        return facets[fi*4 + lv];
    }

    inline int &Quads::vert(const int fi, const int lv) {
        assert(fi>=0 && fi<nfacets() && lv>=0 && lv<4);
        return facets[fi*4 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Polygons::nfacets() const {
        return static_cast<int>(offset.size())-1;
    }

    inline int Polygons::facet_size(const int fi) const {
        assert(fi>=0 && fi<nfacets());
        return offset[fi+1]-offset[fi];
    }

    inline int Polygons::corner(const int fi, const int ci) const {
        assert(fi>=0 && fi<nfacets());
        return offset[fi]+ci;
    }

    inline int Polygons::vert(const int fi, const int lv) const {
        assert(fi>=0 && fi<nfacets());
        assert(lv>=0 && lv<facet_size(fi));
        return facets[offset[fi]+lv];
    }

    inline int &Polygons::vert(const int fi, const int lv) {
        assert(fi>=0 && fi<nfacets());
        assert(lv>=0 && lv<facet_size(fi));
        return facets[offset[fi]+lv];
    }
}

#endif //__SURFACE_H__

