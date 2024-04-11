#ifndef __SPARSEMATRIX_H__
#define __SPARSEMATRIX_H__

#include "linexpr.h"
#include <algorithm>
#include <iostream>
#include <vector>

namespace UM {
    namespace Linear {
        struct RCSMatrix {
            std::vector<LinExpr::Term> mat = {};
            std::vector<int> offset = {0};
        };

        struct LILMatrix {
            std::vector<LinExpr> lines;
        };

        inline LILMatrix lil_identity(int n) {
            LILMatrix m;
            m.lines.reserve(n);
            for (int i=0; i<n; i++)
                m.lines.push_back(X(i));
        }
    }
}

#endif //__SPARSEMATRIX_H__

