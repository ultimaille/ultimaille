#ifndef __PERMUTATION_H__
#define __PERMUTATION_H__

#include <vector>
#include <cassert>

struct Permutation {
    Permutation(const std::vector<int> &p) : destinations(p) {
        is_valid();
    }

    bool is_valid() {
        int n = destinations.size();
        std::vector<bool> visited(n, false);
        for (int i=0; i<n; i++) {
            if (destinations[i] >= n || visited[i])
                return false;
            visited[destinations[i]] = true;
        }
        return true;
    }


    // Permutation method takes A = [a, b, c, d, e] and P = [4, 3, 2, 0, 1], then it returns [e, d, c, a, b].
    // We are allowed to use only constant space and linear time.
    // We can swap each element in A with the right element required by P,
    // after each swap, there will be one more element in the right position,
    // and do this in a circular fashion for each of the positions (swap elements pointed with ^s):
    // [a, b, c, d, e] <- P[0] = 4 != 0 (where a initially was), swap 0 (where a is) with 4
    //  ^           ^
    // [e, b, c, d, a] <- P[4] = 1 != 0 (where a initially was), swap 4 (where a is) with 1
    //     ^        ^
    // [e, a, c, d, b] <- P[1] = 3 != 0 (where a initially was), swap 1 (where a is) with 3
    //     ^     ^
    // [e, d, c, a, b] <- P[3] = 0 == 0 (where a initially was), finish step
    // After one circle, we find the next element in the array that does not stay in the right position, and do this again.
    // So in the end you will get the result you want, and since each position is touched a constant time
    // (for each position, at most one operation (swap) is performed), it is O(n) time.

    // You can store the information of which one is in its right place by:
    // * set the corresponding entry in P to -1, which is unrecoverable: after the operations above,
    //   P will become [-1, -1, 2, -1, -1], which denotes that only the second one might be not in the right position,
    //   and a further step will make sure it is in the right position and terminates the algorithm;
    // * set the corresponding entry in P to -n - 1: P becomes [-5, -4, 2, -1, -2], which can be recovered in O(n) trivially.

    template <typename T> inline void apply(std::vector<T>& data) {
        int n = destinations.size();
        assert(data.size()==static_cast<size_t>(n));
        std::vector<bool> marked(n, false);

        for (int i=0; i<n; ++i) {
//            if (destinations[i] < 0) continue;
            if (marked[i]) continue;
            int cur = i;
            while (destinations[cur] != i) {
//                assert(destinations[cur]>=0);
                assert(!marked[cur]);
                const int target = destinations[cur];
                std::swap(data[cur], data[target]);
                marked[cur] = true;
//                destinations[cur] = -1 - target;
                cur = target;
            }
//           destinations[cur] = -1 - destinations[cur];
            marked[cur] = true;
        }
//      for (int i=0; i<n; ++i)
//          destinations[i] = -1 - destinations[i];
    }

    const std::vector<int> &destinations;
};

#endif //__PERMUTATION_H__

