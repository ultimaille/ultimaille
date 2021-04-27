#ifndef __DISJOINTSET_H__
#define __DISJOINTSET_H__
#include <vector>
#include <numeric>
#include <cassert>

namespace UM {
    // see https://en.wikipedia.org/wiki/Disjoint-set_data_structure#Disjoint-set_forests

    struct DisjointSet {
        // constructor knowing the number of elements
        DisjointSet(const int n) : m_ids(n), m_size(n, 1) {
            std::iota(m_ids.begin(), m_ids.end(), 0); // initialize the union-find data structure (m_ids[i]=i)
        }

        void merge(const int a, const int b) {
            assert(a>=0 && b>=0 && a<size() && b<size());
            int rootA = root(a);
            int rootB = root(b);
            if (rootA == rootB) return;
            if (m_size[rootB] < m_size[rootA]) std::swap(rootA, rootB);
            m_ids[rootA] = rootB;
            m_size[rootB] += m_size[rootA];
        }

        // connect a node and its parents to the root and return the root
        int root(const int i) {
            assert(i>=0 && i<size());
            // find the root
            int root_id = m_ids[i];
            while (root_id!=m_ids[root_id]) root_id = m_ids[root_id];
            // connect all the path to root
            int i_ = i;
            while (i_!=root_id) {
                int new_i = m_ids[i_];
                m_ids[i_] = root_id;
                i_ = new_i;
            }
            return root_id;
        }

        //! returns the size of the set containing i
        int setsize(const int i) {
            assert(i>=0 && i<size());
            return m_size[root(i)];
        }

        //! returns true if the two elements are in the same set
        bool same(const int a, const int b) {
            assert(a>=0 && b>=0 && a<size() && b<size());
            return root(a)==root(b);
        }

        int size() const {
            return m_ids.size();
        }

        // return the number of sets
        int nsets() {
            int num = 0;
            for (int i=0; i<(int)m_ids.size(); i++)
                num += (i == root(i)); // count the different roots
            return num;
        }

        // return the number of set and fill setId
        int get_sets_id(std::vector<int> &id2setid) {
            id2setid.resize(m_ids.size()); // prepare the correspondance root id => set id
            int nsets = 0;
            for (int i=0; i<(int)m_ids.size(); i++) {
                if (i != root(i)) continue;
                id2setid[i] = nsets;
                nsets++;
            }
            for (int i=0; i<(int)m_ids.size(); i++) // fill the rest
                id2setid[i] = id2setid[m_ids[i]];
            return nsets;
        }

        std::vector<int> m_ids;
        std::vector<int> m_size;
    };
}

#endif //__DISJOINTSET_H__

