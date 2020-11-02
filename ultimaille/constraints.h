#ifndef __CONSTRAINTS_H__
#define __CONSTRAINTS_H__
#include <vector>
#include <numeric>
#include <cassert>

namespace UM {
    struct SignedPairwiseEquality {
        SignedPairwiseEquality(const int num) : m_ids(num), m_size(num, 1), m_same_sign(num, true), m_conflicts(num, false) {
            std::iota(m_ids.begin(), m_ids.end(), 0); // initialize the union-find data structure (m_ids[i]=i)
        }

        void merge(const int a, const int b, const bool same_sign) {
            assert(a>=0 && b>=0 && a<size() && b<size());
            int rootA = root(a);
            int rootB = root(b);
            if (rootA == rootB) { // check that the new condition is consistent with the other constraints
                if ((m_same_sign[a]==m_same_sign[b])!=same_sign) // after calling root() a and b are direct children of rootA and rootB
                    m_conflicts[rootA] = true;
                return;
            }

            if (m_size[rootB] < m_size[rootA]) std::swap(rootA, rootB);
            m_ids[rootA] = rootB;

            m_same_sign[rootA] = (m_same_sign[a]==m_same_sign[b])==same_sign;
            m_size[rootB] += m_size[rootA];
            m_conflicts[rootB] = m_conflicts[rootB] | m_conflicts[rootA];
        }

        int root(const int i) { // connect a node and its ancestors to the root and return the root
            assert(i>=0 && i<size());
            bool same_sign_with_root = m_same_sign[i];
            int root_id = i;
            while (root_id!=m_ids[root_id]) {
                same_sign_with_root = (same_sign_with_root==m_same_sign[m_ids[root_id]]);
                root_id = m_ids[root_id];
            }

            // connect all the path to root
            int i_ = i;
            while (i_!=root_id) {
                int new_i = m_ids[i_];
                bool new_same_sign = (same_sign_with_root==m_same_sign[i_]);
                m_ids[i_] = root_id;
                m_same_sign[i_] = same_sign_with_root;
                i_ = new_i;
                same_sign_with_root = new_same_sign;
            }
            return root_id;
        }

        int size() const {
            return m_ids.size();
        }

        int setsize(const int i) {
            assert(i>=0 && i<size());
            return m_size[root(i)];
        }

        int nsets() {
            int num = 0;
            for (int i=0; i<(int)m_ids.size(); i++)
                num += (i == root(i)); // count the different roots
            return num;
        }

        int reduce(std::vector<int> &id, std::vector<int> &sign) {
            int nsets = 0;

            id.resize(size()); // prepare the correspondance root id => set id
            for (int i=0; i<size(); i++) {
                if (i != root(i)) continue;
                id[i] = nsets;
                nsets++;
            }

            for (int i=0; i<(int)m_ids.size(); i++) // fill the rest
                id[i] = id[m_ids[i]];

            sign = std::vector<int>(size(), 0);
            for (int i=0; i<size(); i++) {
                if (m_conflicts[root(i)]) continue;
                sign[i] = (m_same_sign[i] ? 1 : -1);
            }
            return nsets;
        }

        std::vector<int> m_ids;  // disjoint-set forest
        std::vector<int> m_size; // tree size (defined for the tree roots only)
        std::vector<bool> m_same_sign; // signs of the equality constraints, can be incoherent
        std::vector<bool> m_conflicts; // if this boolean is set, the tree signs are incoherent => the variables =0
    };
}

#endif //__CONSTRAINTS_H__

