#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <iostream>
#include <map>

#include "inspect.h"
#include "helpers/disjointset.h"
#include "algebra/covariance.h"
#include "algebra/eigen.h"

namespace UM {
    BBox3 Inspect<PointSet>::bbox() const {
        BBox3 bbox;
        for (vec3 const &p : pts)
            bbox.add(p);
        return bbox;
    }

    vec3 Inspect<PointSet>::barycenter() const {
        vec3 ave = {0, 0, 0};
        for (const vec3 &p : pts)
            ave += p;
        return ave / static_cast<double>(pts.size());
    }

    std::tuple<mat3x3,vec3,vec3> Inspect<PointSet>::principal_axes() const {
        PointSetCovariance cov(*pts.data);
        auto [eval, evec] = eigendecompose_symmetric(cov.cov);
        if (pts.size()<4) return {mat3x3::identity(), {1.,1.,1.}, cov.center}; // If the system is under-determined, return the trivial basis
        return { evec, eval, cov.center };
    }

    double Inspect<Surface>::avg_edge_len() const {
        double sum = 0;
        for (auto f : m.iter_facets())
            for (auto h : f.iter_halfedges())
                sum += vec3(h).norm();
        return sum / static_cast<double>(m.ncorners());
    }

    int Inspect<Surface>::euler_characteristic() const {
        um_assert(m.connected());
        um_assert(is_manifold(false));

        int nb_vertices = 0; // number of **non-isolated** vertices
        for (auto v : m.iter_vertices())
            nb_vertices += v.halfedge().active();

        int nb_facets = 0; // number of *active* facets
        for ([[maybe_unused]] auto f : m.iter_facets())
            nb_facets++;

        int nb_edges = 0;
        for (auto f : m.iter_facets())
            for (auto h : f.iter_halfedges())
                nb_edges += !h.opposite().active() || h.from() > h.to(); // count each edge once
        return nb_vertices + nb_facets - nb_edges;
    }

    int Inspect<Surface>::nb_connected_components() const {
        um_assert(m.connected());
        um_assert(is_manifold(false));

        DisjointSet ds(m.nverts()+1);
        for (auto h : m.iter_halfedges())
            ds.merge(h.from(), h.to());

        for (auto v : m.iter_vertices())
            if (!v.halfedge().active())
                ds.merge(v, m.nverts());
        std::vector<int> id2setid;
        return ds.get_sets_id(id2setid) - 1;
    }

    int Inspect<Surface>::nb_boundaries() const {
        um_assert(m.connected());
        um_assert(is_manifold(false));

        DisjointSet ds(m.nverts()+1);
        for (auto v : m.iter_vertices())
            if (!v.on_boundary() || !v.halfedge().active())
                ds.merge(v, m.nverts());

        for (auto h : m.iter_halfedges()) {
            auto opp = h.opposite();
            if (!opp.active())
                ds.merge(h.from(), h.to());
        }
        std::vector<int> id2setid;
        return ds.get_sets_id(id2setid) - 1;
    }

    bool Inspect<Surface>::is_disk(bool verbose) const {
        um_assert(m.connected());

        if (!is_manifold(verbose))
            return false;

        int nb_comp = nb_connected_components();
        int xi      = euler_characteristic();
        int nb_bnd  = nb_boundaries();

        bool disk = nb_comp==1 && xi==1 && nb_bnd==1;

        if (verbose) {
            std::string message;
            message += "#connected_components = " + std::to_string(nb_comp) + (nb_comp == 1 ? " (expected)" : " (not a disk)") + "\n";
            message += "#boundaries = " + std::to_string(nb_bnd) + (nb_bnd == 1 ? " (expected)" : " (not a disk)") + "\n";
            message += "Euler characteristic = " + std::to_string(xi) + (xi == 1 ? " (expected)" : " (not a disk)") + "\n";
            message += "------------------------------------------------------\n";
            message += std::string("Verdict: the surface is ") + (disk ? "" : "*NOT* ") + "a topological disk\n";
        }

        return disk;
    }

    bool Inspect<Surface>::is_manifold(bool verbose) const {
        um_assert(m.connected());
        int nb_facets = 0; // number of *active* facets
        for ([[maybe_unused]] auto f : m.iter_facets())
            nb_facets++;

        int nb_vertices = m.nverts();

        int nb_isolated_vertices = 0;
        for (auto v : m.iter_vertices())
            nb_isolated_vertices += !v.halfedge().active();

        int nb_null_edges = 0;
        for (auto h : m.iter_halfedges())
            nb_null_edges += h.to() == h.from();

        std::map<int, int> nb_opposites;
        std::map<int, int> nb_duplicates;
        for (auto h : m.iter_halfedges()) {
            int nb_opp = 0;
            int nb_dup = 0;
            for (auto cir : h.from().iter_halfedges()) {
                if (cir.prev().from() == h.to())    nb_opp++;
                if (cir != h && h.to() == cir.to()) nb_dup++;
            }

            if (nb_opposites.contains(nb_opp))
                nb_opposites[nb_opp]++;
            else
                nb_opposites[nb_opp] = 1;

            if (nb_duplicates.contains(nb_dup))
                nb_duplicates[nb_dup]++;
            else
                nb_duplicates[nb_dup] = 1;
        }

        std::vector<int> multiple_umbrella_vertices; // An umbrella of a vertex is a subset of facets incident to the vertex whose adjacency graph is a cycle (for interior vertices, for the boundary vertices it is a path).
        for (auto v : m.iter_vertices()) {
            int nb_halfedges = 0;
            for ([[maybe_unused]] auto h : v.iter_halfedges())
                nb_halfedges++;
            int nb_sector_halfedges = 0;
            if (nb_halfedges)
                for ([[maybe_unused]] auto h : v.halfedge().iter_sector_halfedges())
                    nb_sector_halfedges++;
            if (nb_halfedges != nb_sector_halfedges)
                multiple_umbrella_vertices.push_back(v);
        }

        bool manifold = true;

        if (verbose) {
            std::string message;
            message += "#vertices = " + std::to_string(nb_vertices) + "\n";
            message += "#facets   = " + std::to_string(nb_facets) + "\n";
            message += "#isolated vertices = " + std::to_string(nb_isolated_vertices) + "\n";

            message += "------------------------------------------------------\n";
            for (auto i : nb_opposites) {
                message += "We found " + std::to_string(i.second) + " halfedge(s) with " + std::to_string(i.first) + " opposites" + (i.first == 0 ? " (boundary)" : (i.first == 1 ? " (manifold)" : " (non manifold)")) + "\n";
                manifold = manifold && (i.first < 2);
            }
            message += "------------------------------------------------------\n";
            for (auto i : nb_duplicates) {
                message += "We found " + std::to_string(i.second) + " halfedge(s) duplicated " + std::to_string(i.first) + " times" + (i.first == 0 ? " (expected)" : " (non manifold)") + "\n";
                manifold = manifold && (i.first == 0);
            }
            message += "#null edges = " + std::to_string(nb_null_edges) + (nb_null_edges == 0 ? " (expected)" : " (WARNING)") + "\n";
            manifold = manifold && (nb_null_edges == 0);

            for (int vid : multiple_umbrella_vertices) {
                message += "vertex id = " + std::to_string(vid) + " has mutiple umbrellas (WARNING) \n";
                manifold = false;
            }
            message += "------------------------------------------------------\n";
            message += std::string("Verdict: the surface is ") + (manifold ? "" : "*NOT* ") + "manifold\n";

            std::cerr << message;
        }
        return manifold;
    }
}

