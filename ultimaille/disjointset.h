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

	/*

	struct SignedPairwiseEquality : DisjointSet {
		SignedPairwiseEquality(const int num) : DisjointSet(num), m_same_sign(num, true), m_conflicts(num, false) { }

		void merge(const int a, const int b, const bool same_sign) {
		    assert(a>=0 && b>=0 && a<size() && b<size());
		    int rootA = root(a);
		    int rootB = root(b);
		    if (rootA == rootB) { // check that the new condition is consistent with the other constraints
		        if ((m_same_sign[a]==m_same_sign[b])!=same_sign)
		            m_conflicts[rootA] = true;
		        return;
		    }

		    if (m_size[rootB] < m_size[rootA]) std::swap(rootA, rootB);
		    m_ids[rootA] = rootB;

		    m_same_sign[rootA] = (m_same_sign[a]==m_same_sign[b])==same_sign;
		    m_size[rootB] += m_size[rootA];
		}
		// connect a node and its parents to the root and return the root
		int root(const int i) {
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
		        m_ids[i_] = root_id;
		        bool new_same_sign = (same_sign_with_root==m_same_sign[i_]);
		        m_same_sign[i_] = same_sign_with_root;
		        i_ = new_i;
		        same_sign_with_root = new_same_sign;
		    }
		    return root_id;
		}

		bool is_zero(const int i) {
		    assert(i>=0 && i<size());
		    return m_conflicts[root(i)];
		}

		int reduce(std::vector<int> &id, std::vector<int> &sign) {
		    int nsets = get_sets_id(id);
		    int nvar = m_ids.size();
		    sign = std::vector<int>(nvar, 0);
		    for (int i=0; i<nvar; i++) {
		        if (is_zero(i)) continue;
		        sign[i] = (m_same_sign[i] ? 1 : -1);
	//    bool s;same(i, root(i), s);            sign[i] = (s ? 1 : -1);
		    }
		    return nsets;
		}

		bool same(const int a, const int b) = delete;
		bool same(const int a, const int b, bool &same_sign) {
		    assert(a>=0 && b>=0 && a<size() && b<size());
		    same_sign = (m_same_sign[a]==m_same_sign[b]);
		    return root(a)==root(b);
		}

		std::vector<bool> m_same_sign;
		std::vector<bool> m_conflicts;
	};
	*/
}
#endif //__DISJOINTSET_H__
