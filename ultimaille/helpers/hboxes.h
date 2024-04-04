#ifndef __HBOXES_H__
#define __HBOXES_H__
#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>
#include "ultimaille/algebra/vec.h"

namespace UM {

    template<int n> struct BBox {
        BBox() {
            for (int i = 0; i < n; i++) {
                min[i] =  std::numeric_limits<double>::max();
                max[i] = -std::numeric_limits<double>::max();
            }
        }

        BBox(vec<n> min, vec<n> max) : min(min), max(max) {}

        void add(const BBox<n> &b) {
            if (b.empty()) return;
            add(b.min);
            add(b.max);
        }

        void add(const vec<n> &p) {
            for (int d=0; d<n; d++) {
                min[d] = std::min<double>(min[d], p[d]);
                max[d] = std::max<double>(max[d], p[d]);
            }
        }

        bool empty() const {
            for (int d=0; d<n; d++)
                if (max[d] < min[d])
                    return true;
            return false;
        }

        bool intersect(const BBox<n> &b) const {
            for (int d=0; d<n; d++)
                if (min[d] > b.max[d] || max[d] < b.min[d])
                    return false;
            return true;
        }

        bool contains(const vec<n> &v) const {
            for (int d=0; d<n; d++)
                if (min[d] > v[d] || max[d] < v[d])
                    return false;
            return true;
        }

        void dilate(double eps) {
            vec<n> t;

            for (int d = 0; d < n; d++) {
                t[d] = eps;
            }
            min = min - t;
            max = max + t;
        }

        vec<n> center() const {
            return (min + max)*.5;
        }

        vec<n> size() const {
            return max - min;
        }

        vec<n> min;
        vec<n> max;
    };


    typedef BBox<1> BBox1;
    typedef BBox<2> BBox2;
    typedef BBox<3> BBox3;

    inline unsigned int mylog2(unsigned int x) {
            unsigned int ans = 0 ;
            while (x>>=1) ans++;
            return ans ;
    }

    /**
     * Store bounding boxes as hierarchical tree of boxes
    */
    template<int n> struct HBoxes {

        HBoxes() {
        }


        HBoxes(std::vector<BBox<n>> const &boxes) {
            init(boxes);
        }



        void init(std::vector<BBox<n>> const &boxes) {
            int nboxes = boxes.size();
            std::vector<vec<n>> G(nboxes);
            // Implicit binary tree: 
            // https://opendatastructures.org/ods-cpp/10_1_Implicit_Binary_Tree.html
            // Use indirection map
            tree_pos_to_org.resize(nboxes);
            for (int b=0; b<nboxes; b++) {
                G[b] = boxes[b].center();
                tree_pos_to_org[b] = b;
            }
    // TODO: Maybe not necessary anymore on windows
    #if defined(_OPENMP) && _OPENMP>=200805
    #pragma omp parallel
    #pragma omp single nowait
    #endif
            sort(G, 0, nboxes);

            // Compute the offset that mark the end index of nodes boxes 
            // and the start index of leaves boxes
            offset = static_cast<int>(std::pow(2., 1. + mylog2(nboxes))) - 1;
            // Tree is resized to the number of boxes (leaves) 
            // + number of added boxes (nodes)
            tree.resize(offset + nboxes);

            // Add leaves at the right place using indirection map
            for (int b=0; b<nboxes; b++)
                tree[offset + b] = boxes[tree_pos_to_org[b]];
            // Add nodes at the right place following the pattern of implicit binary tree
            for (int i=offset; i--;) {
                for (int son = 2*i+1; son<2*i+3; son++)
                    if (son < static_cast<int>(tree.size()))
                        tree[i].add(tree[son]);
            }
        }
        // Sort a slice of the indirection array
        void sort(std::vector<vec<n>> &G, int org, int dest) const {
            // Create a bounding box that encompass underlying bounding boxes 
            // (found in the selected slice of indirection array)
            BBox<n> b;
            for (int i=org; i<dest; i++)
                b.add(G[tree_pos_to_org[i]]);

            // Search for the best dim to cut (the largest)
            int dim = n - 1; // find the best dim to cut
            for (int d=0; d<n; d++)
                if (b.max[d] - b.min[d] > b.max[dim] - b.min[dim])
                    dim = d;

            // Sort the slice of the indirection array:
            // order boxes by their center on the split axis
            std::sort(tree_pos_to_org.begin()+org, tree_pos_to_org.begin()+dest, [&G,dim](const int A, const int B){ return G[A][dim]>G[B][dim]; });

            // If it still just one box, we reach the leaf
            if (dest - org <= 2) return;

            // Try to split array (tree) into two parts: 
            // If the tree is not balanced (size : 2^n) we have to find the highest power of 2
            // |---|---|---|
            // 0   16 32   48
            int m = org + static_cast<int>(pow(2., mylog2(dest-org-1)));
    #if defined(_OPENMP) && _OPENMP>=200805
    #pragma omp task
    #endif
            // Sort nodes boxes
            sort(G, org, m);
    #if defined(_OPENMP) && _OPENMP>=200805
    #pragma omp task
    #endif
            // Sort leaves boxes
            sort(G, m, dest);
        }
        
        void intersect(BBox<n> const &b, std::vector<int> &primitives, int node=0) const {
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

        int offset = -1;
        mutable std::vector<int> tree_pos_to_org = {};
        std::vector<BBox<n>> tree = {};
    };

    typedef HBoxes<1> HBoxes1;
    typedef HBoxes<2> HBoxes2;
    typedef HBoxes<3> HBoxes3;
}

#endif //__HBOXES_H__

