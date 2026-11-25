#include <iostream>
#include <algorithm>
#include <cassert>
#include "attributes.h"
#include "polyline.h"
#include "surface.h"
#include "volume.h"

namespace UM {
    int PolyLine::nverts() const {
        return points.size();
    }

    int PolyLine::nedges() const {
        assert(edges.size()%2==0);
        return edges.size()/2;
    }

    int PolyLine::vert(const int s, const int lv) const {
        assert(s>=0 && s<nedges() && lv>=0 && lv<2);
        return edges[s*2 + lv];
    }

    int &PolyLine::vert(const int s, const int lv) {
        assert(s>=0 && s<nedges() && lv>=0 && lv<2);
        return edges[s*2 + lv];
    }

    int PolyLine::create_edges(const int n) {
        edges.resize(edges.size()+n*2);
        resize_attrs();
        return nedges()-n;
    }

    void PolyLine::resize_attrs() {
        for (auto &wp : attr)  if (auto spt = wp.lock())
            spt->resize(nedges());
    }

    void PolyLine::delete_vertices(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nverts());
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

    void PolyLine::delete_edges(const std::vector<bool> &to_kill) {
        assert(!connected());
        assert(to_kill.size()==(size_t)nedges());

        std::vector<int> edges_old2new(nedges(),  -1);
        int new_nb_endpoints = 0;

        for (int s=0; s<nedges(); s++) {
            if (!to_kill[s]) continue;
            edges_old2new[s] = new_nb_endpoints/2;
            for (int lv=0; lv<2; lv++)
                edges[new_nb_endpoints++] = vert(s, lv);
        }

        edges.resize(new_nb_endpoints);

        std::erase_if(attr, [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; });
        for (auto &wp : attr) if (auto spt = wp.lock())
            spt->compress(edges_old2new);
    }

    void PolyLine::connect() {
        if (!conn) conn = std::make_unique<Connectivity>(*this);
        else conn->init();
    }
    void PolyLine::disconnect() {
        conn.reset();
    }

    void PolyLine::compact(bool delete_isolated_vertices) {
        if (!conn) return;
        um_assert(conn->active.ptr != nullptr);
        std::vector<bool> to_kill = conn->active.ptr->data;
        to_kill.flip();
        disconnect(); // happy assert(conn==nullptr)!
        delete_edges(to_kill);
        if (delete_isolated_vertices)
            PolyLine::delete_isolated_vertices();
        connect();
    }

    void PolyLine::delete_isolated_vertices() {
        std::vector<bool> to_kill(nverts(), true);
        for (int e = 0; e < nedges(); e++)
            for (int lv = 0; lv < 2; lv++)
                to_kill[vert(e, lv)] = false;
        delete_vertices(to_kill);
    }

    PolyLine::Connectivity::Connectivity(PolyLine& m) : m(m), v2e(m, -1), e2e(m, -1), active(m, true) {
        init();
    }

    void PolyLine::Connectivity::init() {
        e2e.fill(-1);
        active.fill(true);
        v2e.fill(-1);
        for (int e=0; e<m.nedges(); e++) {
            int org = m.vert(e, 0);
            if (v2e[org] != -1) e2e[e] = v2e[org];
            v2e[org] = e;
        }

        active.fill(true);
        v2e.fill(-1);
        for (int e = 0; e < m.nedges(); e++) {
            int v = m.vert(e, 0);
            e2e[e] = v2e[v];
            v2e[v] = e;
        }
    }

    PolyLine::Edge PolyLine::Connectivity::create_edge(int v0, int v1){
        int e = m.create_edges(1);
        m.vert(e, 0) = v0;
        m.vert(e, 1) = v1;
        um_assert(e2e[e] == -1); //e2e[e] = -1;
        um_assert(active[e]==true);
        int org = m.vert(e, 0);
        e2e[e] = v2e[org];
        v2e[org] = e;
        return Edge(m, e);
    }

    PolyLine::Edge PolyLine::Edge::opposite() {
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

