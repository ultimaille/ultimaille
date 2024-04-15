#ifndef __RBUILDER_H__
#define __RBUILDER_H__

#include "linexpr.h"
#include "sparsematrix.h"

namespace UM {
    namespace Linear {
        struct ReductionBuilder {
            ReductionBuilder(int n) { C = lol_identity(n); }

            void add_constraint(LinExpr &e);
            void add_constraint(LinExpr &&e) { add_constraint(e); }

            // N.B. the builder is destroyed
            CRSMatrix tocrs() {
                CRSMatrix result = C.tocrs();
                C = {};
                return result;
            }

        private:
            // re-express e in terms of free variables
            void reduce(LinExpr &e);

            LOLMatrix C; // square matrix of constraints C x = 0
        };

    }
}

#endif //__RBUILDER_H__

