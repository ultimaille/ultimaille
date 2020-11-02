#ifndef __PERMUTATION_H__
#define __PERMUTATION_H__

#include <vector>
#include <cassert>

namespace UM {
    // TODO: move the permutation data inside the struct (instead of an external ref)
    // TODO: invert the permutation in place (or apply_reverse)
    struct Permutation {
        Permutation(const std::vector<int> &p) : destinations(p) {
            is_valid();
        }

        bool is_valid() {
            int n = destinations.size();
            std::vector<bool> visited(n, false);
            for (int i=0; i<n; i++) {
                if (destinations[i]<0 || destinations[i] >= n || visited[i])
                    return false;
                visited[destinations[i]] = true;
            }
            return true;
        }

        // Imagine a data array A = [a, b, c, d, e] and a permutation P = [4, 3, 2, 0, 1]; after applying P we get [e, d, c, a, b].
        // A trivial solution would be:
        //
        //  template<typename T> void apply_permutation(std::vector<T>& A, const std::vector<int>& P) {
        //      std::vector<T> A2(A.size());
        //      for (size_t i = 0; i < A.size(); i++)
        //          A2[i] = A[P[i]];
        //      A = A2;
        //  }
        //
        // Let us say we want to use constant space and linear time to apply the permutation.
        //
        // P[X] = Y means "an item on position Y should be at position X".
        // So we should move an item that is now at position X somewhere else - swap it with item on position Y.
        // Then we have a right item on position X, but the original X-item now on position Y,
        // maybe should be occupied by someone else (an item Z).
        // So we check P[Y] = Z and move the X-item further until we got P[?] = X which mean that on position ? should be an item
        // from position X - which is exactly the X-item we've been kicking around all this time. Loop closed.

        // Thus, we can swap each element in A with the right element required by P,
        // after each swap, there will be one more element in the right position.
        // We do this in a circular fashion for each of the positions (swap elements pointed with ^s):
        // [a, b, c, d, e] <- P[0] = 4 != 0 (where a initially was), swap 0 (where a is) with 4
        //  ^           ^
        // [e, b, c, d, a] <- P[4] = 1 != 0 (where a initially was), swap 4 (where a is) with 1
        //     ^        ^
        // [e, a, c, d, b] <- P[1] = 3 != 0 (where a initially was), swap 1 (where a is) with 3
        //     ^     ^
        // [e, d, c, a, b] <- P[3] = 0 == 0 (where a initially was), finish step
        // Once this loop terminates, we find the next element in the array that does not stay in the right position, and do this again.

        // The tricky part is to determine whether the element is in its right place.
        // If we wanted to perform the permutation in constant space, we can mark the indices:
        // set the corresponding entry in P to -n - 1: P becomes [-5, -4, 2, -1, -2], which can be recovered in O(n) trivially.
        // It is not thread-safe, thus I prefer to allocate a vector of booleans to indicate the elements in the right place.

#if 0
        // thread-unsafe, but constant memory
        template <typename T> inline void apply(std::vector<T>& data) {
            int n = destinations.size();
            assert(data.size()==static_cast<size_t>(n));
            std::vector<int>& destinations = const_cast<std::vector<int>&>(this->destinations);


            for (int i=0; i<n; ++i) {
                if (destinations[i] < 0) continue;
                int cur = i;
                while (destinations[cur] != i) {
                    assert(destinations[cur]>=0);
                    const int target = destinations[cur];
                    std::swap(data[cur], data[target]);
                    destinations[cur] = -1 - target;
                    cur = target;
                }
                destinations[cur] = -1 - destinations[cur];
            }
            for (int i=0; i<n; ++i)
                destinations[i] = -1 - destinations[i];
        }
#else
        template <typename T> inline void apply(std::vector<T>& data) const {
            int n = destinations.size();
            assert(data.size()==static_cast<size_t>(n));
            std::vector<bool> marked(n, false);

            for (int i=0; i<n; ++i) {
                if (marked[i]) continue;
                int cur = i;
                while (destinations[cur] != i) {
                    assert(!marked[cur]);
                    const int target = destinations[cur];
                    std::swap(data[cur], data[target]);
                    marked[cur] = true;
                    cur = target;
                }
                marked[cur] = true;
            }
        }
#endif

        const std::vector<int> &destinations;
    };
}

#endif //__PERMUTATION_H__

