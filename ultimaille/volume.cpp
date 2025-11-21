#include <cassert>
#include "helpers/disjointset.h"
#include "volume.h"
#include "surface.h"
#include "attr_binding.h"

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
        for (auto cell : m.iter_cells())
            for (auto corner : cell.iter_corners())
                v2c[corner.vertex()] = corner;

        // chain cell corners around vertices
        for (auto cell : m.iter_cells())
            for (auto corner : cell.iter_corners()) {
                c2c[corner] = v2c[corner.vertex()];
                v2c[corner.vertex()] = corner;
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

    Volume::Connectivity::Connectivity(Volume &m) : m(m),  adjacent(m, -1) {
          reset();
    }

    bool are_facets_adjacent(const Volume &m, int c1, int c2, int lf1, int lf2) {
        int n = m.facet_size(lf1);
        if (n!=m.facet_size(lf2)) return false;

        for (int i=0; i<n; i++) {
            bool found = true;
            for (int j=0; found && j<n; j++)
                found = (m.facet_vert(c1, lf1, (i+j)%n) == m.facet_vert(c2, lf2, n-j-1));
            if (found) return true;
        }

        return false;
    }


    void Volume::Connectivity::reset() {
        adjacent.fill(-1); // -1 not found already, -2 more than one candidate

        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(m, v2c, c2c);

        for (int c1=0; c1<m.ncells(); c1++) // TODO rewrite it with iterators
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                int f1 = m.cell(c1).facet(lf1);
                int crn1 = v2c[m.facet_vert(c1, lf1, 0)];
                int crn2 = crn1;
                do {
                    int c2 = crn2 / m.nverts_per_cell();
                    if (c2!=c1) {
                        for (int lf2=0; lf2<m.nfacets_per_cell(); lf2++) {
                            if (!are_facets_adjacent(m, c1, c2, lf1, lf2)) continue;
                            int f2 = m.cell(c2).facet(lf2);
                            if (adjacent[f1]==-1 && adjacent[f2]==-1) { // perfect match
                                adjacent[f1] = f2;
                                adjacent[f2] = f1;
                                continue;
                            }
                            // destroy previous links
                            if (adjacent[f1]>=0) adjacent[adjacent[f1]] = -2;
                            if (adjacent[f2]>=0) adjacent[adjacent[f2]] = -2;

                            adjacent[f1] = -2;
                            adjacent[f2] = -2;
                        }
                    }
                    crn2 = c2c[crn2];
                } while (crn2!=crn1 && adjacent[f1]<0);
            }
        for (int c1=0; c1<m.ncells(); c1++)
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                int f1 = m.cell(c1).facet(lf1);
                if (adjacent[f1] == -2) adjacent[f1] = -1;
            }
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

