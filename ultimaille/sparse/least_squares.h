#ifndef __LEASTSQUARES_H__
#define __LEASTSQUARES_H__
#include <memory>

#include "linexpr.h"
#include "nullspace.h"

namespace UM {

    // least squares solver: basically an OpenNL wrapper
    struct LeastSquares {
        LeastSquares(int nvars, bool verbose = false, double threshold = 1e-10, int nb_max_iter = 5000);
        ~LeastSquares();

        void add_to_energy(const LinExpr& e);
        void fix(int var, double value);
        void solve();
        double value(int i) { return X[i]; }
        int nvars() { return X.size(); }

        std::vector<double> X = {};
    protected:
        void* context = nullptr;
        bool locking_impossible = false; // OpenNL requires to lock all variables before composing the energy
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////

    struct ConstrainedLeastSquares {
        ConstrainedLeastSquares(int nvars, bool verbose = false, double threshold = 1e-10, int nb_max_iter = 5000);
        void add_to_constraints(const LinExpr& le);
        void add_to_energy(const LinExpr& le);
        void solve();

        double value(int i) {
            double dot = 0;
            for (const SparseElement &e : M.iter_row(i))
                dot += e.value * lsptr->X[e.index];
            return dot;
        }

        bool verbose;
        double threshold;
        int nb_max_iter;
        int nfree;
        std::unique_ptr<LeastSquares> lsptr;
        CRSMatrix M;
        NullSpaceBuilder rb;
    };
}

#endif //__LEASTSQUARES_H__

