#ifndef __VOLUME_CONNECTIVITY_H__
#define __VOLUME_CONNECTIVITY_H__

#include <vector>
#include <memory>
#include "syntactic-sugar/assert.h"
#include "polyline.h"

namespace UM {
    struct Volume;

    struct OppositeFacet { // adjacency between half-facets
        OppositeFacet(const Volume &m);
        void reset();

        // [[deprecated]] int & operator[](const int i);
        int   operator[](const int i) const;

        // [[deprecated]] int opposite_c(const int he) const;
        // [[deprecated]] std::vector<int> halfedges_around_edge(const int he) const;

        const Volume &m;
        std::vector<int> adjacent;
    };

}

#endif //__VOLUME_CONNECTIVITY_H__

