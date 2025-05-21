#ifndef __NULLSPACE_H__
#define __NULLSPACE_H__

#include "matrix.h"

namespace UM {
    struct NullSpaceBuilder {

        enum Status {
            SUCCESS    = 0,
            REDUNDANT  = 1,
            IMPOSSIBLE = 2
        };
        
        NullSpaceBuilder(int n, bool free_last=true) : C{lol_identity(n)}, free_last{free_last} { }

        // re-express v in terms of free variables
        // after the call (v.empty()) indicates a redundant constraint;
        // (free_last && !v.empty() && v.front().index == size()-1) indicates an impossible constraint
        void leading_to_free(SparseVector &v);
        Status test(SparseVector &e);
        Status test(SparseVector &&e) { return test(e); }
        void add_constraint(SparseVector &e);
        void add_constraint(SparseVector &&e) { return add_constraint(e); }

        // N.B. the builder is destroyed
        CRSMatrix to_crs() {
            for (SparseVector &r : C.rows)
                leading_to_free(r);
            C.drop_zero_columns();
            CRSMatrix result = C.to_crs();
            C = {};
            return result;
        }

        inline int size() { return C.rows.size(); }

        LOLMatrix C;    // square matrix of constraints C x = 0
        bool free_last; // specify whether last variable can be leading or must be free
    };
}

#endif //__NULLSPACE_H__

