#ifndef __KNN_H__
#define __KNN_H__

#include <utility>
#include <algorithm>
#include <numeric>
#include <queue>
#include <vector>
#include "geometry.h"

namespace UM {
    template<int D> struct KNN { // do not try anything but D=2 or D=3
        const std::vector<vec<D>> &pts;
        const int n;
        std::vector<int> tree;

        KNN(const std::vector<vec<D>> &points): pts(points), n(points.size()), tree(points.size()) {
            std::iota(tree.begin(), tree.end(), 0);
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel
#pragma omp single nowait
#endif
            build(0, n);
        }

        void build(const int L, const int R, const int dim=0) { // build is O(n log n) using divide and conquer
            if (L >= R) return;
            int M = (L+R)/2; // get median in O(n), split dim coordinate
            std::nth_element(tree.begin()+L, tree.begin()+M, tree.begin()+R, [this, dim](int a, int b) { return pts[a][dim]<pts[b][dim]; });
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp task
#endif
            build(L,   M, (dim+1)%D);
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp task
#endif
            build(M+1, R, (dim+1)%D);

        }

        // k-nearest neighbor query, O(k log(k) log(n)) on average
        std::vector<int> query(const vec<D> &p, const int k = 1) const {
            std::priority_queue<std::pair<double, int> > pq; // priority queue for KNN, keep the K nearest

            perform_query(pq, p, k, 0, n, false); // recursion
            std::vector<int> neighbors;
            while (!pq.empty()) { // collect points
                neighbors.push_back(pq.top().second);
                pq.pop();
            }
            std::reverse(neighbors.begin(), neighbors.end());
            return neighbors;
        }

        void perform_query(std::priority_queue<std::pair<double, int> > &pq, const vec<D> &p, const int k, const int L, const int R, const int dim) const {
            if (L >= R) return;
            int M = (L+R)/2;
            vec<D> d = p - pts[tree[M]];
            if (static_cast<int>(pq.size()) < k || d*d < pq.top().first) { // if point is nearer to the kth farthest, put point in queue
                pq.push(std::make_pair(d*d, tree[M]));
                if (static_cast<int>(pq.size()) > k) pq.pop(); // keep k elements only
            }
            int nearL = L, nearR = M, farL = M + 1, farR = R;
            if (d[dim] > 0) { // right is nearer
                std::swap(nearL, farL);
                std::swap(nearR, farR);
            }
            // query the nearer child
            perform_query(pq, p, k, nearL, nearR, (dim+1)%D);

            if (static_cast<int>(pq.size()) < k || d[dim]*d[dim] < pq.top().first) // query the farther child if there might be candidates
                perform_query(pq, p, k, farL, farR, (dim+1)%D);
        }
    };
}

#endif //__KNN_H__

