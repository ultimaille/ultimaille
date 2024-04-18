#ifndef __SPARSE_H__
#define __SPARSE_H__

#include "vector.h"

namespace UM {
    // read-only compressed row storage matrix
    struct CRSMatrix {
        double dot(const std::vector<double> &x) const;
        CRSMatrix transpose() const;

        inline int nrows() const { return offset.size()-1; }
        inline int nnz() const { return offset.back(); }
        int count_columns() const;

        std::vector<SparseElement> mat = {};
        std::vector<int> offset = { 0 };
    };

    // list of lists matrix
    struct LOLMatrix {
        CRSMatrix to_crs();
        void drop_zero_columns();

        inline       SparseVector& operator[](int i)       { return rows[i]; }
        inline const SparseVector& operator[](int i) const { return rows[i]; }
        inline int nrows() const { return rows.size(); }

        // sort and aggregate terms, remove near-zero entries
        void compact();

        // the matrix is supposed to be compacted beforehands
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

