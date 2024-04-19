#ifndef __LINEXPR_H__
#define __LINEXPR_H__

#include "vector.h"

namespace UM {

    // Sparse representation of a linear expression b + \sigma_i a_i x_i.
    // A linear expression consists of a constant term, plus a list of coefficient-variable pairs that capture the linear terms.
    // Linear expressions are used to build linear objective and constraints.
    // They are temporary objects that typically have short lifespans.

    struct LinExpr : SparseVector {
        LinExpr() : SparseVector() {}
        LinExpr(SparseElement e) : SparseVector(e) {}
        LinExpr(std::initializer_list<SparseElement> l) : SparseVector(l) {}
        LinExpr(double v) : SparseVector({-1, v}) {}
        LinExpr(SparseVector &&v) { data = std::move(v.data); } // N.B. no compact() call here
    };

    inline LinExpr operator-(const LinExpr& a, const LinExpr& b) {
        return (SparseVector)a - (SparseVector)b;
    }

    inline LinExpr operator+(const LinExpr& a, const LinExpr& b) {
        return (SparseVector)a + (SparseVector)b;
    }

    inline LinExpr operator*(const LinExpr& v, double a) {
        return (SparseVector)v * a;
    }

    inline LinExpr operator*(double a, const LinExpr& v) {
        return v * a;
    }

    namespace Linear {
        inline LinExpr X(int i) { return {{i, 1.}}; }
    }
}

#endif //__LINEXPR_H__

