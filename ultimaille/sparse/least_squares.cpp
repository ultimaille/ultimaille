#include <OpenNL_psm/OpenNL_psm.h>
#include "least_squares.h"
#include "../syntactic-sugar/assert.h"
#include <iostream>

namespace UM {
    LeastSquares::LeastSquares(int nvars, bool verbose, double threshold, int nb_max_iter) : X(nvars, 0.) {
        context = nlNewContext();
        if (verbose)
            nlEnable(NL_VERBOSE);
        nlEnable(NL_VARIABLES_BUFFER);
        nlSolverParameterd(NL_THRESHOLD, threshold);
        nlSolverParameteri(NL_MAX_ITERATIONS, nb_max_iter);
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlSolverParameteri(NL_NB_VARIABLES, (NLint)nvars);
        nlBegin(NL_SYSTEM);
        nlBindBuffer(NL_VARIABLES_BUFFER, 0, (void*)X.data(), (NLuint)sizeof(double));
    }

    LeastSquares::~LeastSquares() {
        if (context) {
            nlDeleteContext(context);
            context = nullptr;
        }
    }

    void LeastSquares::fix(int var, double value) {
        um_assert(!locking_impossible);
        X[var] = value;
        nlMakeCurrent(context);
        nlLockVariable(NLint(var));
    }

    void LeastSquares::add_to_energy(const LinExpr& le) {
        nlMakeCurrent(context);
        if (!locking_impossible) {
            nlBegin(NL_MATRIX);
            locking_impossible = true;
        }
        nlBegin(NL_ROW);
        double rhs = 0.;
        for (const SparseElement& e : le)
            if (e.index<0) rhs -= e.value;
            else nlCoefficient(e.index, e.value);
        nlRightHandSide(rhs);
        nlEnd(NL_ROW);
    }

    void LeastSquares::solve() {
        nlMakeCurrent(context);
        if (!locking_impossible)
            nlBegin(NL_MATRIX);
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();
        nlDeleteContext(context);
        context = nullptr;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////

    ConstrainedLeastSquares::ConstrainedLeastSquares(int nvars, bool verbose, double threshold, int nb_max_iter) :
        verbose(verbose), threshold(threshold), nb_max_iter(nb_max_iter), nfree(-1), lsptr{nullptr}, M{}, rb(nvars+1, true) {
    }

    static SparseVector le2sp(const LinExpr& le, int const_id) {
        SparseVector v = le;
        v.compact();
        if (!v.empty() && v.front().index<0) {
            v.front().index = const_id;
            v.compact();
        }
        return v;
    }

    void ConstrainedLeastSquares::add_to_constraints(const LinExpr& le) {
        um_assert(!lsptr);
        rb.add_constraint(le2sp(le, rb.C.nrows()-1));
    }

    void ConstrainedLeastSquares::add_to_energy(const LinExpr& le) {
        if (!lsptr) {
            M = rb.to_crs();
            nfree = M.count_columns()-1;
            lsptr = std::make_unique<LeastSquares>(nfree+1, verbose, threshold, nb_max_iter);
            lsptr->fix(nfree, 1.);
        }
        SparseVector v = le;
        v.compact();
        if (!v.empty() && v.front().index<0)
            v.front().index = M.nrows()-1;
        lsptr->add_to_energy(v * M);
    }

    void ConstrainedLeastSquares::solve() {
        lsptr->solve();
    }

}

