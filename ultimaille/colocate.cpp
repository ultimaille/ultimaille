#include "knn.h"
#include "colocate.h"

namespace UM {
    void colocate(const std::vector<vec3> &points, std::vector<int> &old2new, double tolerance) {
        int nb = points.size();
        old2new = std::vector<int>(nb, -1);

        KNN<3> knn(points);

#pragma omp parallel for
        for (int seed=0; seed<nb; seed++) {
            int k = std::min<int>(6, nb);
            while (1) {
                std::vector<int> neighbors = knn.query(points[seed], k);

                bool allfound = false;
                int smallest = seed;
                for (int i=0; i<k; i++) {
                    double dist2 = (points[seed]-points[neighbors[i]]).norm2();
                    if (dist2 > tolerance*tolerance) {
                        allfound = true;
                        break;
                    }
                    smallest = std::min<int>(smallest, neighbors[i]);
                }
                old2new[seed] = smallest;
                if (allfound || k==nb) break;

                k = std::min<int>(k+k/2, nb);
            }
        }

        for (int i=0; i<nb; i++) {
            int j = i;
            // colocate clusters of identical vertices onto smallest index
            while (old2new[j] != j)
                j = old2new[j];

            old2new[i] = j;
        }
    }
}

