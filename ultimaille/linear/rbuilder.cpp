#include "rbuilder.h"
#include "../syntactic-sugar/assert.h"

namespace UM {
    namespace Linear {
        void ReductionBuilder::reduce(LinExpr &e) {
            LinExpr result;
            for (const LinExpr::Term &t : e) {                    // for each term in the linear expression
                LinExpr &row = C[t.i];
                if (t.i<0 || (row.size()==1 && row[0].i==t.i))    // if constant or free variable,
                    result.data.push_back(t);                     // push directly into the resulting expression
                else {                                            // otherwise (if basic variable)
                    reduce(row);                                  // reduce the constraint itself
                    for (const LinExpr::Term &s : row)            // and push the scaled constraint into the result
                        result.data.emplace_back(s.i, t.a * s.a);
                }
            }
            result.compact();                                     // do not forget to aggregate terms
            e = std::move(result);
        }

        void ReductionBuilder::add_constraint(LinExpr &e) {
            reduce(e);
            if (e.empty()) return;
            um_assert(e.back().i>=0); // check for the impossible constraint non-zero constant = 0
            LinExpr::Term last = e.back();
            e.data.pop_back();
            e *= -1. / last.a;
            C[last.i] = std::move(e);
        }
    }
}

