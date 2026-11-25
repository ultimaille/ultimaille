#ifndef __VOLUME_H__
#define __VOLUME_H__
#include <vector>
#include <array>
#include <memory>
#include <cassert>

#include "attributes.h"
#include "pointset.h"
#include "volume_reference.h"
#include "polyline.h"

namespace UM {
    struct Volume {
        enum CELL_TYPE { TETRAHEDRON=0, HEXAHEDRON=1, WEDGE=2, PYRAMID=3 };
        CELL_TYPE cell_type;

        PointSet points = {};
        std::vector<int> cells = {};

        std::vector<std::weak_ptr<ContainerBase>> attr_cells   = {};
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

        template <typename T> void delete_vertices(const T&to_kill);
        template <typename T> void delete_cells(const T &to_kill);
        void delete_isolated_vertices();
        int  create_cells(const int n);

        int nverts()     const;
        int ncells()     const;
        int nfacets()    const;
        int ncorners()   const;
        int nhalfedges() const;

        constexpr int     nverts_per_cell() const;
        constexpr int    nfacets_per_cell() const;
        constexpr int nhalfedges_per_cell() const;
        constexpr int facet_size(const int f) const;

        int facet_vert(const int c, const int lf, const int lv) const;
        int  vert(const int c, const int lv) const;
        int &vert(const int c, const int lv);

        Volume(CELL_TYPE cell_type) : cell_type(cell_type) {}
        Volume(const Volume& m)            = delete;
        Volume(Volume&& m)                 = delete;
        Volume& operator=(const Volume& m) = delete;

        template <typename T> struct CellAttribute;
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
        struct Corner;
        struct Halfedge;
        struct Facet;
        struct Cell;

        Vertex     vertex(int id) const;
        Corner     corner(int id) const;
        Halfedge halfedge(int id) const;
        Facet       facet(int id) const;
        Cell         cell(int id) const;

        auto iter_vertices()  const;
        auto iter_corners()   const;
        auto iter_halfedges() const;
        auto iter_facets()    const;
        auto iter_cells()     const;

        inline bool connected() const;
        void connect();
        void disconnect();

        std::unique_ptr<Connectivity> conn = {};
        private:
        void compute_corner_to_corner_map(std::vector<int> &v2c, std::vector<int> &c2c) const;
        void resize_attrs();
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Tetrahedra : Volume {
        Tetrahedra() : Volume(Volume::TETRAHEDRON) {}
    };

    struct Hexahedra : Volume {
        Hexahedra() : Volume(Volume::HEXAHEDRON) {}
    };

    struct Wedges : Volume {
        Wedges() : Volume(Volume::WEDGE) {}
    };

    struct Pyramids : Volume {
        Pyramids() : Volume(Volume::PYRAMID) {}
    };

    /*
     * EdgeGraph represents edges of a volumetric mesh
     *    -> Concept: An edge corresponds to a set of halfedges with same from/to vertices
     *    -> Design choice: Edges are manifold, so a non manifold edge is represented edges with same from/to vertices
     *   -> Benefits: A) use EdgeAttributes and B) gives acces to vertex neigborhood i.e. all primitives can reach its neigborgs
     */

    struct EdgeGraph : public PolyLine {
        EdgeGraph(Volume& m);
        PolyLine::Edge edge_from_halfedge(Volume::Halfedge h);
        Volume::Halfedge halfedge_from_edge(Edge e);

        std::vector<int> m_halfedge_from_edge;
        Volume& m;
    };


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // these implementations are here and not in the .cpp because all inline functions must be available in all translation units //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline int Volume::nverts() const {
        return points.size();
    }

    inline int Volume::ncorners() const {
        return cells.size();
    }

    inline int Volume::ncells() const {
        return cells.size() / nverts_per_cell();
    }

    inline int Volume::nfacets() const {
        return ncells() * nfacets_per_cell();
    }

    inline int Volume::nhalfedges() const {
        return ncells() * nhalfedges_per_cell();
    }

    inline int Volume::vert(const int c, const int lv) const {
        assert(c>=0 && c<ncells() && lv>=0 && lv<nverts_per_cell());
        return cells[c*nverts_per_cell() + lv];
    }

    inline int &Volume::vert(const int c, const int lv) {
        assert(c>=0 && c<ncells() && lv>=0 && lv<nverts_per_cell());
        return cells[c*nverts_per_cell() + lv];
    }

    inline constexpr int Volume::nverts_per_cell() const {
        return reference_cells[cell_type].nverts();
    }

    inline constexpr int Volume::nfacets_per_cell() const {
        return reference_cells[cell_type].nfacets();
    }

    inline constexpr int Volume::nhalfedges_per_cell() const {
        return reference_cells[cell_type].ncorners();
    }

    inline constexpr int Volume::facet_size(const int f) const {
        assert(f>=0 && f<nfacets());
        return reference_cells[cell_type].facet_size(f % nfacets_per_cell());
    }

    inline int Volume::facet_vert(const int c, const int lf, const int lv) const {
        assert(c>=0 && c<ncells() && lf>=0 && lf<nfacets_per_cell() && lv>=0 && lv<facet_size(lf));
        return vert(c, reference_cells[cell_type].vert(lf, lv));
    }

#include "volume_attributes.h"
#include "volume_primitives.h"
#include "volume_iterators.h"

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    template <typename T> void UM::Volume::delete_cells(const T &to_kill) {
        assert(!connected());

        constexpr bool invocable = std::is_invocable_r_v<bool, T, int>;
        if constexpr (!invocable)
            assert(to_kill.size()==(size_t)ncells());

        std::vector<int>   cells_old2new(ncells(),   -1);
        std::vector<int>  facets_old2new(nfacets(),  -1);
        std::vector<int> corners_old2new(ncorners(), -1);

        int new_nb_cells   = 0;
        int new_nb_facets  = 0;
        int new_nb_corners = 0;

        for (auto c : iter_cells()) {
            if constexpr (invocable) {
                if (to_kill(c)) continue;
            } else {
                if (to_kill[c]) continue;
            }

            for (auto f : c.iter_facets())
                facets_old2new[f] = new_nb_facets++;
            for (auto corner : c.iter_corners()) {
                corners_old2new[corner] = new_nb_corners;
                cells[new_nb_corners] = corner.vertex();
                new_nb_corners++;
            }
            cells_old2new[c] = new_nb_cells++;
        }

        cells.resize(new_nb_corners);

        std::erase_if(attr_cells,   [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; }); // remove dead attributes
        std::erase_if(attr_facets,  [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; });
        std::erase_if(attr_corners, [](std::weak_ptr<ContainerBase> ptr) { return ptr.lock()==nullptr; });

        for (auto &wp : attr_cells)   if (auto spt = wp.lock())
            spt->compress(cells_old2new);
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->compress(facets_old2new);
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->compress(corners_old2new);
    }

    template <typename T> void Volume::delete_vertices(const T &to_kill) {
        constexpr bool invocable = std::is_invocable_r_v<bool, T, int>;
        if constexpr (!invocable)
            assert(to_kill.size() == (size_t)nverts());

        std::vector<bool> cells_to_kill(ncells(), false);
        std::vector<int> c2c, v2c;
        compute_corner_to_corner_map(v2c, c2c);

        for (int v=0; v<nverts(); v++) {
            if constexpr (invocable) {
                if (!to_kill(v)) continue;
            } else {
                if (!to_kill[v]) continue;
            }
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
}

#endif //__VOLUME_H__

