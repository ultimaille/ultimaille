#ifndef __SPARSE_H__
#define __SPARSE_H__

#include "vector.h"

namespace UM {
    // read-only compressed row storage matrix
    struct CRSMatrix {
        // compute Ax+b
        double dot(const std::vector<double> &x) const;
        //          void transpose();

        inline int nrows() const { return offset.size()-1; }
        std::vector<SparseElement> mat = {};
        std::vector<int> offset = { 0 };
    };

    // list of lists matrix
    struct LOLMatrix {
        CRSMatrix tocrs();
        void drop_zero_columns();

        inline       SparseVector& operator[](int i)       { return rows[i]; }
        inline const SparseVector& operator[](int i) const { return rows[i]; }
        inline int nrows() const { return rows.size(); }

        // sort and aggregate terms, remove near-zero entries
        void compact();

        int count_columns() const;
        int count_nnz() const;

        std::vector<SparseVector> rows = {};
    };

    inline LOLMatrix lol_identity(int n) {
        LOLMatrix m;
        m.rows.reserve(n);
        for (int i=0; i<n; i++)
            m.rows.push_back({{1, 1.}});
        return m;
    }
}

#endif //__SPARSE_H__

