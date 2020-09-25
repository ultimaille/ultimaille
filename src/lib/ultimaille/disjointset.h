#ifndef __DISJOINTSET_H__
#define __DISJOINTSET_H__
#include <vector>
#include <numeric>

// see https://en.wikipedia.org/wiki/Disjoint-set_data_structure#Disjoint-set_forests

struct DisjointSet {
    // constructor knowing the number of elements
    DisjointSet(const int num) : m_ids(num), m_size(num, 1) {
        std::iota(m_ids.begin(), m_ids.end(), 0); // initialize the union-find data structure (m_ids[i]=i)
    }

    void merge(const int a, const int b) {
        int rootA = root(a);
        int rootB = root(b);
        if (rootA == rootB) return;
        if (m_size[rootB] < m_size[rootA]) std::swap(rootA, rootB);
        m_ids[rootA] = rootB;
        m_size[rootB] += m_size[rootA];
    }

    // connect a node and its parents to the root and return the root
    int root(const int i) {
        int i_ = i;
        // find the root
        int rootId = m_ids[i_];
        while (rootId!=m_ids[rootId]) rootId = m_ids[rootId];
        // connect all the path to root
        while (i_!=rootId) {
            int newI = m_ids[i_];
            m_ids[i_] = rootId;
            i_ = newI;
        }
        return rootId;
    }

    //! returns the size of the set containing i
    int setsize(const int i) {
        return m_size[root(i)];
    }

    //! returns true if the two elements are in the same set
    bool same(const int a, const int b) {
        return root(a)==root(b);
    }

    // return the number of sets
    int nsets() {
        int num = 0;
        for (int i=0; i<(int)m_ids.size(); i++) {
            if (i == root(i)) ++num; // count the different roots
        }
        return num;
    }

    // return the number of set and fill setId
    int get_sets_id(std::vector<int> &id2setid) {
        // prepare the correspondance root Id => setId
        id2setid.resize(m_ids.size());
        int nsets = 0;
        for (int i=0; i<(int)m_ids.size(); i++) {
            if (i == root(i)) {
                id2setid[i] = nsets;
                nsets++;
            }
        }
        // fill the other set Ids
        for (int i=0; i<(int)m_ids.size(); i++)
            id2setid[i] = id2setid[m_ids[i]];
        return nsets;
    }

    std::vector<int> m_ids;
    std::vector<int> m_size;
};

#endif //__DISJOINTSET_H__

