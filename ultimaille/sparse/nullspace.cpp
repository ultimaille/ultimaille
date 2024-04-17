#include "nullspace.h"
#include "../syntactic-sugar/assert.h"

namespace UM {
    bool NullSpaceBuilder::reduce(SparseVector &v) {
        SparseVector result;
        for (const SparseElement &e : v) {                             // for each term of the constraint
            SparseVector &row = C[e.index];
            if (row.size()==1 && row[0].index==e.index)                // if variable is free,
                result.data.push_back(e);                              // push the term directly into the resulting expression
            else {                                                     // otherwise
                reduce(row);                                           // reduce the basic variable constraint itself
                for (const SparseElement &e2 : row)                    // and push the scaled constraint into the result
                    result.data.emplace_back(e.index, e.value * e2.value);
            }
        }
        result.compact();                                              // do not forget to aggregate terms
        v = std::move(result);
        return !free_last || v.empty() || v.front().index < size()-1;  // check whether it would be ok to add the constraint
    }

    void NullSpaceBuilder::add_constraint(SparseVector &e) {
        um_assert(reduce(e));                                          // check for the impossible constraint non-zero constant = 0
        if (e.empty()) return;
        SparseElement pivot = free_last && e.back().index==size()-1 ?  // if we need to keep the last variable free,
                              e[rand() % (e.size()-1)] :               // it cannot be used as a pivot
                              e[rand() % e.size()];                    // aside this, the choice of pivot is arbitrary
        e -= SparseVector{pivot};
        e *= -1. / pivot.value;
        C[pivot.index] = std::move(e);
    }
}

