#ifndef __SURFACE_H__
#define __SURFACE_H__
#include <vector>
#include <memory>
#include "attributes.h"
#include "pointset.h"

namespace UM {

    struct Surface { // polygonal mesh interface
        PointSet points = {};
        std::vector<int> facets = {};
        std::vector<std::weak_ptr<ContainerBase>> attr_facets  = {};
        std::vector<std::weak_ptr<ContainerBase>> attr_corners = {};

        ////////////////////////////////////////////////////
        //       _                 _                _     //
        //      | | _____      __ | | _____   _____| |    //
        //      | |/ _ \ \ /\ / / | |/ _ \ \ / / _ \ |    //
        //      | | (_) \ V  V /  | |  __/\ V /  __/ |    //
        //      |_|\___/ \_/\_/   |_|\___| \_/ \___|_|    //
        //                                                //
        ////////////////////////////////////////////////////

        void delete_vertices(const std::vector<bool>& to_kill);
        virtual void delete_facets(const std::vector<bool>& to_kill);
        void delete_isolated_vertices();
        void resize_attrs(); // TODO hide it
        void compress_attrs(const std::vector<bool>& facets_to_kill); // TODO remove this

        int nverts()   const;
        int ncorners() const;

        virtual int nfacets() const = 0;
        virtual int facet_size(const int fi) const = 0;
        virtual int corner(const int fi, const int ci) const = 0;
        virtual int   vert(const int fi, const int lv) const = 0;
        virtual int&  vert(const int fi, const int lv) = 0;

        Surface() = default;
        Surface(const Surface&  m) = delete;
        Surface(      Surface&& m) = delete;
        Surface& operator=(const Surface& m) = delete;

        template <typename T> struct FacetAttribute;
        template <typename T> struct CornerAttribute;

        //////////////////////////////////////////////////////////
        //       _     _       _       _                _       //
        //      | |   (_)     | |     | |              | |      //
        //      | |__  _  __ _| |__   | | _____   _____| |      //
        //      | '_ \| |/ _` | '_ \  | |/ _ \ \ / / _ \ |      //
        //      | | | | | (_| | | | | | |  __/\ V /  __/ |      //
        //      |_| |_|_|\__, |_| |_| |_|\___| \_/ \___|_|      //
        //                __/ |                                 //
        //               |___/                                  //
        //////////////////////////////////////////////////////////

        struct Connectivity;
        struct Primitive;
        struct Vertex;
        struct Halfedge;
        struct Facet;

        Vertex     vertex(int id);
        Halfedge halfedge(int id);
        Facet       facet(int id);

        auto iter_vertices();
        auto iter_halfedges();
        auto iter_facets();

        inline bool connected() const;
        void connect();
        void disconnect();
        void compact(bool delete_isolated_vertices = true);

        std::unique_ptr<Connectivity> conn = {};
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Triangles : Surface { // simplicial mesh implementation
        Triangles() = default;
        int create_facets(const int n);

        int nfacets() const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int   vert(const int fi, const int lv) const;
        int&  vert(const int fi, const int lv);
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Quads : Surface { // quad mesh implementation
        Quads() = default;
        int create_facets(const int n);

        int nfacets()  const;
        int facet_size(const int) const;
        int corner(const int fi, const int ci) const;
        int   vert(const int fi, const int lv) const;
        int&  vert(const int fi, const int lv);
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Polygons : Surface { // polygonal mesh implementation
        std::vector<int> offset = { 0 };
        Polygons() = default;

        int create_facets(const int n, const int size);
        void delete_facets(const std::vector<bool>& to_kill);

        int nfacets()  const;
        int facet_size(const int fi) const;
        int corner(const int fi, const int ci) const;
        int   vert(const int fi, const int lv) const;
        int&  vert(const int fi, const int lv);
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Surface::nverts() const {
        return points.size();
    }

    inline int Surface::ncorners() const {
        return facets.size();
    }

    inline bool Surface::connected() const {
        return conn != nullptr;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Triangles::nfacets() const {
        assert(0 == facets.size() % 3);
        return facets.size() / 3;
    }

    inline int Triangles::facet_size(const int) const {
        return 3;
    }

    inline int Triangles::corner(const int fi, const int ci) const {
        assert(ci >= 0 && ci < 3 && fi >= 0 && fi < nfacets());
        return fi * 3 + ci;
    }

    inline int Triangles::vert(const int fi, const int lv) const {
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 3);
        return facets[fi * 3 + lv];
    }

    inline int& Triangles::vert(const int fi, const int lv) {
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 3);
        return facets[fi * 3 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Quads::nfacets() const {
        assert(0 == facets.size() % 4);
        return facets.size() / 4;
    }

    inline int Quads::facet_size(const int) const {
        return 4;
    }

    inline int Quads::corner(const int fi, const int ci) const {
        assert(ci >= 0 && ci < 4 && fi >= 0 && fi < nfacets());
        return fi * 4 + ci;
    }

    inline int Quads::vert(const int fi, const int lv) const {
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 4);
        return facets[fi * 4 + lv];
    }

    inline int& Quads::vert(const int fi, const int lv) {
        assert(fi >= 0 && fi < nfacets() && lv >= 0 && lv < 4);
        return facets[fi * 4 + lv];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Polygons::nfacets() const {
        return static_cast<int>(offset.size()) - 1;
    }

    inline int Polygons::facet_size(const int fi) const {
        assert(fi >= 0 && fi < nfacets());
        return offset[fi + 1] - offset[fi];
    }

    inline int Polygons::corner(const int fi, const int ci) const {
        assert(fi >= 0 && fi < nfacets());
        return offset[fi] + ci;
    }

    inline int Polygons::vert(const int fi, const int lv) const {
        assert(fi >= 0 && fi < nfacets());
        assert(lv >= 0 && lv < facet_size(fi));
        return facets[offset[fi] + lv];
    }

    inline int& Polygons::vert(const int fi, const int lv) {
        assert(fi >= 0 && fi < nfacets());
        assert(lv >= 0 && lv < facet_size(fi));
        return facets[offset[fi] + lv];
    }

}

#include "surface_attributes.h"
#include "surface_primitives.h"
#include "surface_iterators.h"

#endif //__SURFACE_H__

