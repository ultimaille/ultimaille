#include <iostream>
#include <algorithm>
#include <cassert>

#include "volume.h"
#include "volume_connectivity.h"
#include "attributes.h"

namespace UM {
    // global halfedge id from local id
    int HalfEdgeHelper::halfedge(const int cell, const int cell_facet, const int facet_he) const {
        assert(cell>=0 && cell<m.ncells());
        assert(cell_facet>=0 && cell_facet<m.nfacets_per_cell());
        assert(facet_he>=0 && facet_he<=m.facet_size(cell_facet));
        return cell*nhalfedges_per_cell() + reference_cells[m.cell_type].corner(cell_facet, facet_he);
    }

    // c, org and dst are global indices
    int HalfEdgeHelper::halfedge_from_verts(const int c, const int org, const int dst) const {
        assert(c>=0 && c<m.ncells());
        assert(org>=0 && org<m.nverts());
        assert(dst>=0 && dst<m.nverts());

        for (int cf=0; cf<m.nfacets_per_cell(); cf++) {
            int nbv = m.facet_size(cf);
            for (int cfh=0; cfh<nbv; cfh++) {
                if (
                        org == m.facet_vert(c, cf, cfh) &&
                        dst == m.facet_vert(c, cf, (cfh+1) % nbv)
                   )
                    return halfedge(c, cf, cfh);
            }
        }
        return -1;
    }

    int HalfEdgeHelper::nhalfedges_per_cell() const {
        return reference_conn[m.cell_type].m.ncorners();
    }

    int HalfEdgeHelper::nhalfedges() const {
        return m.ncells() * nhalfedges_per_cell();
    }

    vec3 HalfEdgeHelper::geom(const int he) const {
        return m.points[to(he)] - m.points[from(he)];
    }

    // global cell id
    int HalfEdgeHelper::cell(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return he / nhalfedges_per_cell();
    }

    // global facet id
    int HalfEdgeHelper::facet(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return m.facet(cell(he), cell_facet(he));
    }

    // global corner id
    int HalfEdgeHelper::corner(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return cell(he) * m.nverts_per_cell() + reference_conn[m.cell_type].from(cell_halfedge(he));
    }

    // local facet id
    int HalfEdgeHelper::cell_facet(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return reference_conn[m.cell_type].facet(cell_halfedge(he));
    }

    // local halfedge id
    int HalfEdgeHelper::cell_halfedge(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return he % nhalfedges_per_cell();
    }

    // local halfedge id
    int HalfEdgeHelper::facet_halfedge(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return cell_halfedge(he) - reference_conn[m.cell_type].m.corner(cell_facet(he), 0);
    }

    // global vertex id
    int HalfEdgeHelper::from(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return m.facet_vert(cell(he), cell_facet(he), facet_halfedge(he));
    }

    // global vertex id
    int HalfEdgeHelper::to(const int he) const {
        assert(he>=0 && he<nhalfedges());
        const int size = m.facet_size(cell_facet(he));
        return m.facet_vert(cell(he), cell_facet(he), (facet_halfedge(he)+1) % size);
    }

    // global cell halfedge id
    int HalfEdgeHelper::prev(const int he) const {
        assert(he>=0 && he<nhalfedges());
        if (facet_halfedge(he)>0) return he - 1;
        const int size = m.facet_size(cell_facet(he));
        return he + size - 1;
    }

    // global cell halfedge id
    int HalfEdgeHelper::next(const int he) const {
        assert(he>=0 && he<nhalfedges());
        const int size = m.facet_size(cell_facet(he));
        if (facet_halfedge(he)<size-1) return he + 1;
        return he - size + 1;
    }

    int HalfEdgeHelper::opposite_f(const int he) const {
        assert(he>=0 && he<nhalfedges());
        return nhalfedges_per_cell()*cell(he) + reference_conn[m.cell_type].opposite(cell_halfedge(he));
    }

    int HalfEdgeHelper::opposite_c(const OppositeFacet &adj, const int he) const {
        return adj.opposite_c(he);
    }

    void compute_corner_to_corner_map(const Volume &m, std::vector<int> &v2c, std::vector<int> &c2c);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    OppositeFacet::OppositeFacet(const Volume &m) : m(m), adjacent(m.nfacets(), -1) {
        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(m, v2c, c2c);

        for (int c1=0; c1<m.ncells(); c1++)
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                int f1 = m.facet(c1, lf1);
                if (adjacent[f1]>=0) continue;
                int crn1 = v2c[m.facet_vert(c1, lf1, 0)];
                int crn2 = crn1;
                do {
                    int c2 = crn2 / m.nverts_per_cell();
                    if (c2!=c1) {
                        for (int lf2=0; lf2<m.nfacets_per_cell(); lf2++) {
                            if (!are_facets_adjacent(m, c1, c2, lf1, lf2)) continue;
                            int f2 = m.facet(c2, lf2);
                            adjacent[f1] = f2;
                            adjacent[f2] = f1;
                        }
                    }
                    crn2 = c2c[crn2];
                } while (crn2!=crn1 && adjacent[f1]<0);
            }
    }

    int & OppositeFacet::operator[](const int i)       { assert(i>=0 && i<m.nfacets()); return adjacent[i]; }
    int   OppositeFacet::operator[](const int i) const { assert(i>=0 && i<m.nfacets()); return adjacent[i]; }

    int OppositeFacet::opposite_c(const int he) const {
        assert(he>=0 && he<m.heh.nhalfedges());
        int hfacet = m.heh.facet(he);
        assert((int)adjacent.size()>hfacet);
        int opp_hfacet = adjacent[hfacet];
        if (opp_hfacet<0) return -1;
        int opp_cell = m.cell_from_facet(opp_hfacet);
        return m.heh.halfedge_from_verts(opp_cell, m.heh.to(he), m.heh.from(he));
    }

    std::vector<int> OppositeFacet::halfedges_around_edge(const int he) const {
        std::vector<int> result;
        int around_e_cir = he;
        do { // rewind if boundary
            result.push_back(around_e_cir);
            if (opposite_c(around_e_cir)<0) break;
            around_e_cir = m.heh.opposite_f(opposite_c(around_e_cir));
        } while (around_e_cir != he);

        if (opposite_c(around_e_cir)>=0 && around_e_cir == he)
            return result;

        result.clear(); // iterate forward if the edge is on border
        do {
            result.push_back(around_e_cir);
            around_e_cir = opposite_c(m.heh.opposite_f(around_e_cir));
        } while (around_e_cir>=0);
        return result;
    }

    halfedge_around_edge_iter::halfedge_around_edge_iter(const OppositeFacet &of, const int he) : of(of), heh(of.m.heh) {
        assert(he>=0 && he<heh.nhalfedges());
        int cur = he;
        do {
            int oppc = of.opposite_c(cur);
            if (oppc<0) {
                start = cur;
                return;
            }
            cur = heh.opposite_f(oppc);
        } while (cur != he);
        start = finish = he;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED  //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VolumeConnectivity::VolumeConnectivity(const Volume &p_m) : m(p_m), adjacent(m.nfacets(), -1) {
        std::cerr << "WARNING: attention, this interface is now deprecated and will be soon removed" << std::endl;
        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(m, v2c, c2c);

        for (int c1=0; c1<m.ncells(); c1++)
            for (int lf1=0; lf1<m.nfacets_per_cell(); lf1++) {
                int f1 = m.facet(c1, lf1);
                if (adjacent[f1]>=0) continue;
                int crn1 = v2c[m.facet_vert(c1, lf1, 0)];
                int crn2 = crn1;
                do {
                    int c2 = crn2 / m.nverts_per_cell();
                    if (c2!=c1) {
                        for (int lf2=0; lf2<m.nfacets_per_cell(); lf2++) {
                            if (!are_facets_adjacent(m, c1, c2, lf1, lf2)) continue;
                            int f2 = m.facet(c2, lf2);
                            adjacent[f1] = f2;
                            adjacent[f2] = f1;
                        }
                    }
                    crn2 = c2c[crn2];
                } while (crn2!=crn1 && adjacent[f1]<0);
            }
    }

    // c, org and dst are global indices
    int VolumeConnectivity::halfedge_from_verts(const int c, const int org, const int dst) const {
        assert(c>=0); assert(org>=0); assert(dst>=0);

        for (int cf=0; cf<m.nfacets_per_cell(); cf++) {
            int nbv = m.facet_size(cf);
            for (int cfh=0; cfh<nbv; cfh++) {
                if (
                        org == m.facet_vert(c, cf, cfh) &&
                        dst == m.facet_vert(c, cf, (cfh+1) % nbv)
                   )
                    return halfedge(c, cf, cfh);
            }
        }
        return -1;
    }

    int VolumeConnectivity::opposite_f(const int he) const {
        assert(he>=0);
        return halfedge_from_verts(cell(he), to(he), from(he));
    }

    int VolumeConnectivity::opposite_c(const int he) const {
        assert(he>=0);
        int hfacet = facet(he);
        assert((int)adjacent.size()>hfacet);
        int opp_hfacet = adjacent[hfacet];
        if (opp_hfacet<0) return -1;
        int opp_cell = m.cell_from_facet(opp_hfacet);
        return halfedge_from_verts(opp_cell, to(he), from(he));
    }

    std::vector<int> VolumeConnectivity::halfedges_around_edge(const int he) const {
        std::vector<int> result;
        int around_e_cir = he;
        do { // rewind if boundary
            result.push_back(around_e_cir);
            if (opposite_c(around_e_cir)<0) break;
            around_e_cir = opposite_f(opposite_c(around_e_cir));
        } while (around_e_cir != he);

        if (opposite_c(around_e_cir)>=0 && around_e_cir == he)
            return result;

        result.clear(); // iterate forward if the edge is on border
        do {
            result.push_back(around_e_cir);
            around_e_cir = opposite_c(opposite_f(around_e_cir));
        } while (around_e_cir>=0);
        return result;
    }

    // global halfedge id from local id
    int VolumeConnectivity::halfedge(const int cell, const int cell_facet, const int facet_he) const {
        assert(cell>=0 && cell<m.ncells());
        assert(cell_facet>=0 && cell_facet<m.nfacets_per_cell());
        assert(facet_he>=0 && facet_he<=m.facet_size(cell_facet));
        return max_f*max_h*cell + max_h*cell_facet + facet_he;
    }

    // global cell id
    int VolumeConnectivity::cell(const int he) const {
        assert(he>=0);
        return he / (max_f*max_h);
    }

    // global facet id
    int VolumeConnectivity::facet(const int he) const {
        assert(he>=0);
        return m.facet(cell(he), cell_facet(he));
    }

    // local facet id
    int VolumeConnectivity::cell_facet(const int he) const {
        assert(he>=0);
        return (he % (max_f*max_h)) / max_h;
    }

    // local halfedge id
    int VolumeConnectivity::facet_halfedge(const int he) const {
        assert(he>=0);
        return he % max_h;
    }

    int VolumeConnectivity::facet_size(const int he) const {
        assert(he>=0);
        return m.facet_size(cell_facet(he));
    }

    // global vertex id
    int VolumeConnectivity::from(const int he) const {
        assert(he>=0);
        return m.facet_vert(cell(he), cell_facet(he), facet_halfedge(he));
    }

    // global vertex id
    int VolumeConnectivity::to(const int he) const {
        assert(he>=0);
        return m.facet_vert(cell(he), cell_facet(he), (facet_halfedge(he)+1)%facet_size(he));
    }

    // global cell halfedge id
    int VolumeConnectivity::prev(const int he) const {
        assert(he>=0);
        if (facet_halfedge(he)>0) return he - 1;
        return he + facet_size(he) - 1;
    }

    // global cell halfedge id
    int VolumeConnectivity::next(const int he) const {
        assert(he>=0);
        const int size = facet_size(he);
        if (facet_halfedge(he)<size-1) return he + 1;
        return he - size + 1;
    }
}

