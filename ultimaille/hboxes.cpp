#include <algorithm>
#include <cassert>
#include <cmath>
#include "hboxes.h"

namespace UM {
    inline unsigned int mylog2(unsigned int x) {
        unsigned int ans = 0 ;
        while (x>>=1) ans++;
        return ans ;
    }

    HBoxes::HBoxes(std::vector<BBox3> const &boxes) {
        int nboxes = boxes.size();
        std::vector<vec3> G(nboxes);
        tree_pos_to_org.resize(nboxes);
        for (int b=0; b<nboxes; b++) {
            G[b] = boxes[b].center();
            tree_pos_to_org[b] = b;
        }

#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel
#pragma omp single nowait
#endif
        sort(G, 0, nboxes);

        offset = std::pow(2., 1. + mylog2(nboxes)) - 1;
        tree.resize(offset + nboxes);

        for (int b=0; b<nboxes; b++)
            tree[offset + b] = boxes[tree_pos_to_org[b]];
        for (int i=offset; i--;) {
            for (int son = 2*i+1; son<2*i+3; son++)
                if (son < static_cast<int>(tree.size()))
                    tree[i].add(tree[son]);
        }
    }

    void HBoxes::sort(std::vector<vec3> &G, int org, int dest) const {
        BBox3 b;
        for (int i=org; i<dest; i++)
            b.add(G[tree_pos_to_org[i]]);

        int dim = 2; // find the best dim to cut
        for (int d=0; d<2; d++)
            if (b.max[d] - b.min[d] > b.max[dim] - b.min[dim])
                dim = d;

        std::sort(tree_pos_to_org.begin()+org, tree_pos_to_org.begin()+dest, [&G,dim](const int A, const int B){ return G[A][dim]>G[B][dim]; });

        if (dest - org <= 2) return;

        int m = org + pow(2., mylog2(dest-org-1));
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp task
#endif
        sort(G, org, m);
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp task
#endif
        sort(G, m, dest);
    }

    void HBoxes::intersect(BBox3 const &b, std::vector<int> &primitives, int node) const {
        if (!node) primitives.resize(0);
        assert(node>=0 && node < static_cast<int>(tree.size()));
        if (!tree[node].intersect(b)) return;
        if (node >= offset)
            primitives.push_back(tree_pos_to_org[node - offset]);
        else {
            for (int son=2*node+1; son<2*node+3; son++)
                if (son < static_cast<int>(tree.size()))
                    intersect(b, primitives, son);
        }
    }
}

