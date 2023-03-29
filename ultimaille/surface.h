#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <vector>
#include <memory>
#include "algebra/vec.h"
#include "attributes.h"
#include "pointset.h"
#include "syntactic-sugar/assert.h"

namespace UM {
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
            disconnect();
        }

        Surface() = default;
        Surface(const Surface& m) {
            um_assert(!m.points.size() && !m.facets.size());
        }
        Surface& operator=(const Surface& m) {
            clear();
            um_assert(!m.points.size() && !m.facets.size());
            return *this;
        }

        struct Util {
            const Surface &m;
            vec3 bary_verts(const int f) const;
        } util = {*this};

        struct Connectivity {
            Connectivity(Surface &mesh);
            PointAttribute<int>  v2c; // vertex to corner
            CornerAttribute<int> c2f; // corner to facet
            CornerAttribute<int> c2c; // corner to next corner sharing the same vertex (unordered)
            FacetAttribute<bool> selection; // facets to keep after compacting
        };

        std::unique_ptr<Connectivity> conn = {};
        inline bool connected() { return conn!=nullptr; }

        void connect();
        void disconnect();
        void compact(bool delete_isolated_vertices = true);

/*
        struct Vertex;
        struct Halfedge;
        struct Facet;

        struct Primitive {
            Primitive(Surface &m, int id);
            Primitive(Primitive &p) = default;

            Primitive& operator=(Primitive &p);
            Primitive& operator=(int i);

            operator int  () const;
            operator int& ();

            protected:

            const Surface &m;
            int id;
        };
*/
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Triangles : Surface { // simplicial mesh implementation
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        Triangles() = default;
        Triangles(const Triangles& m) = default;
        Triangles& operator=(const Triangles& m) {
            Surface::operator=(m);
            return *this;
        }

        struct Util : Surface::Util {
            double unsigned_area(const int f) const;
            void project(const int t, vec2& z0, vec2& z1, vec2& z2) const;
            vec3 normal(const int f) const;
        } util = {*this};
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Quads : Surface { // quad mesh implementation
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int  vert(const int fi, const int lv) const;
        int &vert(const int fi, const int lv);

        Quads() = default;
        Quads(const Quads& m) = default;
        Quads& operator=(const Quads& m) {
            Surface::operator=(m);
            return *this;
        }

        struct Util : Surface::Util {
            double unsigned_area(const int f) const;
            vec3 normal(const int f) const;
        } util = {*this};
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Polygons : Surface { // polygonal mesh implementation
        std::vector<int> offset = {0};
        Polygons() = default;
        Polygons(const Polygons& m) = default;
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
        } util = {*this};
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

