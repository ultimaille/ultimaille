#include "nullspace.h"
#include "../syntactic-sugar/assert.h"

namespace UM {
    void NullSpaceBuilder::reduce(SparseVector &v) {
        SparseVector result;
        for (const SparseElement &e : v) {                             // for each term of the constraint
            SparseVector &row = C[e.index];
            if (row.size()==1 && row[0].index==e.index)                // if variable is free,
                result.data.push_back(e);                              // push the term directly into the resulting expression
            else {                                                     // otherwise
                reduce(row);                                           // reduce the basic variable constraint itself
                for (const SparseElement &e2 : row)                    // and push the scaled constraint into the result
                    result.data.push_back(e2 * e.value);
            }
        }
        result.compact();                                              // do not forget to aggregate terms
        v = std::move(result);
    }

    void NullSpaceBuilder::add_constraint(SparseVector &v) {
        reduce(v);

        if (v.empty()) return;
        um_assert(!free_last || v.front().index < size()-1);           // check for the impossible constraint non-zero constant = 0
        SparseElement pivot = free_last && v.back().index==size()-1 ?  // if we need to keep the last variable free,
                              v[rand() % (v.size()-1)] :               // it cannot be used as a pivot
                              v[rand() % v.size()];                    // aside this, the choice of pivot is arbitrary
        v -= SparseVector{pivot};
        v *= -1. / pivot.value;
        C[pivot.index] = std::move(v);
    }
}

