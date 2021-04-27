#ifndef __PERMUTATION_H__
#define __PERMUTATION_H__

#include <vector>
#include <cassert>

namespace UM {
    struct Permutation {
        Permutation(int n) : ind(n) {
            std::iota(ind.begin(), ind.end(), 0);
        }

        int size() const { return ind.size(); }
        auto begin() { return ind.begin(); }
        auto end()   { return ind.end();   }
        int &operator[](const int i) { return ind[i]; }

        bool is_valid() {
            int n = size();
            std::vector<bool> visited(n, false);
            for (int i=0; i<n; i++) {
                if (ind[i]<0 || ind[i] >= n || visited[ind[i]])
                    return false;
                visited[ind[i]] = true;
            }
            return true;
        }

        // Imagine a data array A = [a, b, c, d, e] and a permutation P = [4, 3, 2, 0, 1]; after applying P we get [e, d, c, a, b].
        // A trivial solution is given in apply_lin function that creates a copy of the data to permute.
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


//      // thread-unsafe, but constant memory
//      template <typename T> inline void apply(std::vector<T>& data) {
//          int n = ind.size();
//          assert(data.size()==static_cast<size_t>(n));
//          std::vector<int>& ind = const_cast<std::vector<int>&>(this->ind);
//          for (int i=0; i<n; ++i) {
//              if (ind[i] < 0) continue;
//              int cur = i;
//              while (ind[cur] != i) {
//                  assert(ind[cur]>=0);
//                  const int target = ind[cur];
//                  std::swap(data[cur], data[target]);
//                  ind[cur] = -1 - target;
//                  cur = target;
//              }
//              ind[cur] = -1 - ind[cur];
//          }
//          for (int i=0; i<n; ++i)
//              ind[i] = -1 - ind[i];
//      }

        template <typename T> void apply(std::vector<T>& data) const {
            int n = size();
            assert(static_cast<int>(data.size())==n);
            std::vector<bool> marked(n, false);
            for (int i=0; i<n; ++i) {
                if (marked[i]) continue;
                int cur = i;
                while (ind[cur] != i) {
                    assert(!marked[cur]);
                    std::swap(data[cur], data[ind[cur]]);
                    marked[cur] = true;
                    cur = ind[cur];
                }
                marked[cur] = true;
            }
        }

        template <typename T> void apply_reverse(std::vector<T>& data) const {
            int n = size();
            assert(static_cast<int>(data.size())==n);
            std::vector<bool> marked(n, false);
            for (int i=0; i<n; ++i) {
                if (marked[i]) continue;
                int cur = i;
                do {
                    assert(!marked[ind[cur]]);
                    std::swap(data[i], data[ind[cur]]);
                    marked[ind[cur]] = true;
                    cur = ind[cur];
                } while (ind[cur] != i);
                marked[i] = true;
            }
        }

        template <typename T> void apply_lin(std::vector<T>& data) const {
            int n = size();
            assert(static_cast<int>(data.size())==n);
            std::vector<T> data2(n);
            for (int i = 0; i<n; i++)
                data2[i] = data[ind[i]];
            data = data2;
        }

        template <typename T> void apply_reverse_lin(std::vector<T>& data) const {
            int n = size();
            assert(static_cast<int>(data.size())==n);
            std::vector<T> data2(n);
            for (int i = 0; i<n; i++)
                data2[ind[i]] = data[i];
            data = data2;
        }

        std::vector<int> ind;
    };
}

#endif //__PERMUTATION_H__

