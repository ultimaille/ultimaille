#ifndef __NULLSPACE_H__
#define __NULLSPACE_H__

#include "matrix.h"

namespace UM {
    struct NullSpaceBuilder {
        NullSpaceBuilder(int n, bool free_last=true) : C{lol_identity(n)}, free_last{free_last} { }

        // re-express v in terms of free variables
        // return value is a check whether the constraint is admissible (last variable can be free if needeed)
        bool reduce(SparseVector &v);

        void add_constraint(SparseVector &e);
        void add_constraint(SparseVector &&e) { add_constraint(e); }

        // N.B. the builder is destroyed
        CRSMatrix tocrs() {
            C.drop_zero_columns();
            CRSMatrix result = C.tocrs();
            C = {};
            return result;
        }

        inline int size() { return C.rows.size(); }

        LOLMatrix C;    // square matrix of constraints C x = 0
        bool free_last; // specify whether last variable can be leading or must be free
    };
}

#endif //__NULLSPACE_H__

