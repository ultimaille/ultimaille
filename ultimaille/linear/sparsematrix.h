#ifndef __SPARSEMATRIX_H__
#define __SPARSEMATRIX_H__

#include "linexpr.h"
#include <vector>

namespace UM {
    namespace Linear {

        // read-only compressed row storage matrix
        struct CRSMatrix {
            // compute Ax+b
            double dot(const std::vector<double> &x) const;

            inline int nrows() const { return offset.size()-1; }
            std::vector<LinExpr::Term> mat = {};
            std::vector<int> offset = { 0 };
        };

        // list of lists matrix
        struct LOLMatrix {
            CRSMatrix tocrs();
            void drop_zero_columns();

            // sort and aggregate terms, remove near-zero entries
            void compact();

            // N.B. counts over the base part of the augmented matrix, i.e. it ignores the affine constant
            int count_base_columns() const;

            // count all non-zero entries, including affine constants
            int count_nnz() const;
            inline int nrows() const { return rows.size(); }

            std::vector<LinExpr> rows = {};
        };

        inline LOLMatrix lol_identity(int n) {
            LOLMatrix m;
            m.rows.reserve(n);
            for (int i=0; i<n; i++)
                m.rows.push_back(X(i));
            return m;
        }
    }
}

#endif //__SPARSEMATRIX_H__

