#ifndef __POLYLINE_H__
#define __POLYLINE_H__

#include <vector>
#include <memory>
#include <type_traits>
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

        template <typename T> void delete_vertices(const T &to_kill);
        template <typename T> void delete_edges(const T &to_kill);

        void delete_isolated_vertices();

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
        private:
        void resize_attrs();
    };
}

#include "polyline_attributes.h"
#include "polyline_primitives.h"
#include "polyline_iterators.h"

namespace UM {
    template <typename T> void PolyLine::delete_edges(const T &to_kill) {
        assert(!connected());
        constexpr bool invocable = std::is_invocable_r_v<bool, T, int>;
        if constexpr (!invocable)
            assert(to_kill.size() == (size_t)nedges());

        std::vector<int> edges_old2new(nedges(),  -1);
        int new_nb_endpoints = 0;
        for (int s=0; s<nedges(); s++) {
            if constexpr (invocable) {
                if (to_kill(s)) continue;
            } else {
                if (to_kill[s]) continue;
            }

            edges_old2new[s] = new_nb_endpoints/2;
            for (int lv=0; lv<2; lv++)
                edges[new_nb_endpoints++] = vert(s, lv);
        }
        edges.resize(new_nb_endpoints);

        std::erase_if(attr, [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; });
        for (auto &wp : attr) if (auto spt = wp.lock())
            spt->compress(edges_old2new);
    }

    template <typename T> void PolyLine::delete_vertices(const T &to_kill) {
        assert(!connected());
        constexpr bool invocable = std::is_invocable_r_v<bool, T, int>;
        if constexpr (!invocable)
            assert(to_kill.size() == (size_t)nverts());

        std::vector<bool> edges_to_kill(nedges(), false);
        for (int s=0; s<nedges(); s++)
            for (int lv=0; lv<2; lv++)
                if (to_kill[vert(s, lv)])
                    edges_to_kill[s] = true;
        delete_edges(edges_to_kill);

        std::vector<int> old2new;
        points.delete_points(to_kill, old2new);
        for (int &v : edges)
            v = old2new[v];
    }
}

#endif //__POLYLINE_H__

