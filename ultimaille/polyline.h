#ifndef __POLYLINE_H__
#define __POLYLINE_H__

#include <vector>
#include <memory>
#include "pointset.h"
#include "attributes.h"
#include "syntactic-sugar/assert.h"

namespace UM {
    struct PolyLine {
        PointSet points = {};
        std::vector<int> edges = {};
        std::vector<std::weak_ptr<ContainerBase>> attr = {};

        ////////////////////////////////////////////////////
        //       _                 _                _     //
        //      | | _____      __ | | _____   _____| |    //
        //      | |/ _ \ \ /\ / / | |/ _ \ \ / / _ \ |    //
        //      | | (_) \ V  V /  | |  __/\ V /  __/ |    //
        //      |_|\___/ \_/\_/   |_|\___| \_/ \___|_|    //
        //                                                //
        ////////////////////////////////////////////////////

        void delete_vertices(const std::vector<bool> &to_kill);
        void delete_edges(const std::vector<bool> &to_kill);
        void delete_isolated_vertices();
        void resize_attrs();
        void compress_attrs(const std::vector<bool> &edges_to_kill);

        int nverts() const;
        int nedges() const;

        int  vert(const int s, const int lv) const;
        int &vert(const int s, const int lv)      ;

        int create_edges(const int n);

        PolyLine() = default;
        PolyLine(const PolyLine&  m) = delete;
        PolyLine(      PolyLine&& m) = delete;
        PolyLine& operator=(const PolyLine& m) = delete;

        template <typename T> struct Attribute;


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
        struct Edge;

        Vertex vertex(int id);
        Edge     edge(int id);

        auto iter_vertices();
        auto iter_edges();

        inline bool connected() const;
        void connect();
        void disconnect();
        void compact(bool delete_isolated_vertices = true);

        std::unique_ptr<Connectivity> conn = {};
    };
}

#include "polyline_attributes.h"                         
#include "polyline_primitives.h"
#include "polyline_iterators.h"         


namespace UM {
    inline PolyLine::Edge PolyLine::Edge::opposite() {
        assert(m.connected());
        assert(active());
        Edge result = { m, -1 }; // not found
        for (auto candidate : to().iter_edges()) {
            um_assert(candidate.from() == to());
            if (candidate.to() == from()) {
                if (result == -1) result = candidate;
                else return { m, -1 }; // found more than one
            }
        }
        return result;
    }
}


#endif //__POLYLINE_H__

