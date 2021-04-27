#ifndef __VOLUME_CONNECTIVITY_H__
#define __VOLUME_CONNECTIVITY_H__
#include <vector>
#include <memory>
#include "syntactic-sugar/assert.h"
#include "volume.h"

namespace UM {
    struct VolumeConnectivity { // half-edge-like connectivity interface
        VolumeConnectivity(const Volume &p_m);

        int halfedge(const int cell, const int cell_facet, const int facet_he) const;
        int halfedge_from_verts(const int c, const int org, const int dst) const;

        int           cell(const int he) const;
        int          facet(const int he) const; // TODO: facet -> halffacet?
        int     cell_facet(const int he) const;
        int facet_halfedge(const int he) const; // TODO find better name
        int     facet_size(const int he) const;
        int           from(const int he) const;
        int             to(const int he) const;
        int           prev(const int he) const;
        int           next(const int he) const;

        int opposite_c(const int he) const; // TODO speed-up the implementation
        int opposite_f(const int he) const; // of these two functions

        std::vector<int> halfedges_around_edge(const int he) const;

        static constexpr int max_f = 6; // max number of faces per cell
        static constexpr int max_h = 4; // max number of halfedges per face
        const Volume &m;
        std::vector<int> adjacent; // adjacency between half-facets
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // these implementations are here and not in the .cpp because all inline functions must be available in all translation units //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // global halfedge id from local id
    inline int VolumeConnectivity::halfedge(const int cell, const int cell_facet, const int facet_he) const {
        assert(cell>=0 && cell<m.ncells());
        assert(cell_facet>=0 && cell_facet<m.nfacets_per_cell());
        assert(facet_he>=0 && facet_he<=m.facet_size(cell, cell_facet));
        return max_f*max_h*cell + max_h*cell_facet + facet_he;
    }

    // global cell id
    inline int VolumeConnectivity::cell(const int he) const {
        assert(he>=0);
        return he / (max_f*max_h);
    }

    // global facet id
    inline int VolumeConnectivity::facet(const int he) const {
        assert(he>=0);
        return m.facet(cell(he), cell_facet(he));
    }

    // local facet id
    inline int VolumeConnectivity::cell_facet(const int he) const {
        assert(he>=0);
        return (he % (max_f*max_h)) / max_h;
    }

    // local halfedge id
    inline int VolumeConnectivity::facet_halfedge(const int he) const {
        assert(he>=0);
        return he % max_h;
    }

    inline int VolumeConnectivity::facet_size(const int he) const {
        assert(he>=0);
        return m.facet_size(cell(he), cell_facet(he));
    }

    // global vertex id
    inline int VolumeConnectivity::from(const int he) const {
        assert(he>=0);
        return m.facet_vert(cell(he), cell_facet(he), facet_halfedge(he));
    }

    // global vertex id
    inline int VolumeConnectivity::to(const int he) const {
        assert(he>=0);
        return m.facet_vert(cell(he), cell_facet(he), (facet_halfedge(he)+1)%facet_size(he));
    }

    // global cell halfedge id
    inline int VolumeConnectivity::prev(const int he) const {
        assert(he>=0);
        if (facet_halfedge(he)>0) return he - 1;
        return he + facet_size(he) - 1;
    }

    // global cell halfedge id
    inline int VolumeConnectivity::next(const int he) const {
        assert(he>=0);
        const int size = facet_size(he);
        if (facet_halfedge(he)<size-1) return he + 1;
        return he - size + 1;
    }
}

#endif //__VOLUME_CONNECTIVITY_H__

