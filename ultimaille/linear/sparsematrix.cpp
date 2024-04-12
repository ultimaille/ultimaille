#include "sparsematrix.h"
#include "linexpr.h"
#include <vector>

namespace UM {
    namespace Linear {

        double CRSMatrix::dot(const std::vector<double> &x) const {
            double result = 0;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for
#endif
            for (int i = 0; i < nrows(); ++i) {
                for (int j = offset[i-1]; j<offset[i]; ++j)
                    if (mat[j].i>=0)
                        result += mat[j].a * x[mat[j].i];
                    else
                        result += mat[j].a;
            }
            return result;
        }

        void LOLMatrix::compact() {
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for
#endif
            for (int i = 0; i < (int)rows.size(); ++i)
                rows[i].compact();

        }

        void LOLMatrix::drop_zero_columns() {
            compact();

            // mark non-zero columns
            int ncols = count_base_columns();
            std::vector<bool> zero_cols(ncols, true);

            for (int i = 0; i < (int)rows.size(); ++i) // vector<bool> isn't thread safe
                for (const LinExpr::Term &t : rows[i].data)
                    if (t.i>=0) // ignore the affine constant
                        zero_cols[t.i] = false;

            // remap the columns
            std::vector<int> old2new(ncols);
            int new_col_count = 0;
            for (int col = 0; col < ncols; ++col)
                if (!zero_cols[col])
                    old2new[col] = new_col_count++;

            // apply the map
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for
#endif
            for (int i = 0; i < (int)rows.size(); ++i)
                for (LinExpr::Term &t : rows[i].data)
                    if (t.i>=0)
                        t.i = old2new[t.i];
        }

        int LOLMatrix::count_base_columns() const {
            int max_index = -1;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(max:max_index)
#endif
            for (int i = 0; i < (int)rows.size(); ++i)
                for (const LinExpr::Term &c : rows[i].data)
                    max_index = std::max(max_index, c.i);
            return max_index + 1;
        }

        int LOLMatrix::count_nnz() const {
            int cnt = 0;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(+:cnt)
#endif
            for (int i = 0; i < (int)rows.size(); ++i)
                cnt += rows[i].data.size();
            return cnt;
        }

        CRSMatrix LOLMatrix::tocrs() {
            compact();
            int nnz = count_nnz();
            CRSMatrix m = {
                std::vector(nnz, LinExpr::Term{}),
                std::vector<int>(nrows(), 0)
            };
            for (int i=1; i<=nrows(); ++i) {
                m.offset[i] = m.offset[i-1];
                for (const LinExpr::Term &c : rows[i].data)
                    m.mat[m.offset[i]++] = c;
            }
            return m;
        }
    }
}

