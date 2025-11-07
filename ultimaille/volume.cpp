#include <cassert>
#include "helpers/disjointset.h"
#include "volume.h"
#include "attributes.h"

namespace UM {

    void compute_corner_to_corner_map(const Volume &m, std::vector<int> &v2c, std::vector<int> &c2c);

    void Volume::resize_attrs() {
        for (auto &wp : attr_cells)   if (auto spt = wp.lock())
            spt->resize(ncells());
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->resize(nfacets());
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->resize(ncorners());
    }

    void Volume::delete_vertices(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nverts());
        std::vector<bool> cells_to_kill(ncells(), false);

        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(*this, v2c, c2c);

        for (int v=0; v<nverts(); v++) {
            if (!to_kill[v]) continue;
            int cir = v2c[v];
            if (cir<0) continue; // isolated vertex
            do {
                cells_to_kill[cir/nverts_per_cell()] = true;
                cir = c2c[cir];
            } while (cir != v2c[v]);
        }
        delete_cells(cells_to_kill);

        std::vector<int> old2new;
        points.delete_points(to_kill, old2new);
        for (int &v : cells)
            v = old2new[v];
    }

    void Volume::delete_isolated_vertices()  {
        std::vector<bool> to_kill(nverts(), true);
        for (int c=0; c<ncells(); c++)
            for (int lv=0; lv<nverts_per_cell(); lv++)
                to_kill[vert(c, lv)] = false;
        delete_vertices(to_kill);
    }

    int Volume::create_cells(const int n) {
        assert(n>=0);
        cells.resize(cells.size()+n*nverts_per_cell());
        resize_attrs();
        return ncells()-n;
    }

    void compute_corner_to_corner_map(const Volume &m, std::vector<int> &v2c, std::vector<int> &c2c) {
        c2c = std::vector<int>(m.ncorners(), -1);
        v2c = std::vector<int>(m.nverts(),   -1);

        // compute vertex-to-cell-corner map
        for (int c=0; c<m.ncells(); c++)
            for (int lv=0; lv<m.nverts_per_cell(); lv++)
                v2c[m.vert(c, lv)] = m.corner(c, lv);

        // chain cell corners around vertices
        for (int c=0; c<m.ncells(); c++)
            for (int lv=0; lv<m.nverts_per_cell(); lv++) {
                int crn = m.corner(c, lv);
                int v = m.vert(c, lv);
                c2c[crn] = v2c[v];
                v2c[v] = crn;
            }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    void Volume::connect() {
        if (!conn) conn = std::make_unique<Connectivity>(*this);
        else conn->reset();
    }

    void Volume::disconnect() {
        conn.reset();
    }

    Volume::Connectivity::Connectivity(Volume &m) : m(m),  oppf(m) {
//        reset();
    }

    void Volume::Connectivity::reset() {
        oppf.reset();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    EdgeGraph::EdgeGraph(Volume& m) : m(m) {
        points = m.points;

        DisjointSet ds(m.nhalfedges());
        for(auto h:m.iter_halfedges())
            if (h.facet().opposite().active())
                ds.merge(h,h.opposite_c().opposite_f());

        std::vector<int> edges;
        int nedges = ds.get_sets_id(edges);
        create_edges(nedges);
        m_halfedge_from_edge.resize(nedges,-1);
        for(auto h:m.iter_halfedges()) {
            int e = edges[h];
            vert(e, 0) = h.from();
            vert(e, 1) = h.to();
            m_halfedge_from_edge[e] = h;
        }
        connect();
    }

    PolyLine::Edge EdgeGraph::edge_from_halfedge(Volume::Halfedge h) {
        Vertex from = vertex(h.from());
        Vertex to   = vertex(h.to());
        PolyLine::Edge res(*this,-1);
        for (auto e_manifold : from.iter_edges()){
            if (e_manifold.to()!=to)
                continue;
            if (res==-1){
                res = e_manifold;
                continue;
            }
            // the edge is non manifold, find the one reachable from h
            for (int cir : h.iter_CCW_around_edge())
                for (auto e : from.iter_edges())
                    if (halfedge_from_edge(e)==cir)
                        return e;
        }
        return res;
    }


}

