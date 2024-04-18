#include "matrix.h"

namespace UM {

    double CRSMatrix::dot(const std::vector<double> &x) const {
        double result = 0;
#pragma omp parallel for reduction(+:result)
        for (int i = 0; i < nrows(); ++i) {
            for (int j = offset[i-1]; j<offset[i]; ++j)
                if (mat[j].index>=0)
                    result += mat[j].value * x[mat[j].index];
                else
                    result += mat[j].value;
        }
        return result;
    }

    int CRSMatrix::count_columns() const {
        int max_index = -1;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(max:max_index)
#endif
        for (int i = 0; i < nnz(); ++i)
            max_index = std::max(max_index, mat[i].index);
        return max_index + 1;
    }

    CRSMatrix CRSMatrix::transpose() const {
        int ncols = count_columns();

        // transposed matrix
        std::vector newmat(nnz(), SparseElement(0, 0.));
        std::vector newoffset(ncols+1, 0);

        // count nnz per column
        for (int row = 0; row < nrows(); row++)
            for (int j = offset[row]; j<offset[row+1]; ++j)
                newoffset[1+mat[j].index]++;

        // sum the counts to obtain the column offsets
        for (int col = 0; col < ncols; col++)
            newoffset[col+1] += newoffset[col];

        // fill in the transposed data
        std::vector<int> cur_nnz_per_col(ncols, 0);
        for (int row = 0; row < nrows(); row++)
            for (int j = offset[row]; j<offset[row+1]; ++j)
                newmat[ newoffset[mat[j].index] + cur_nnz_per_col[mat[j].index]++ ] = { row, mat[j].value };

        return { newmat, newoffset };
    }

    void LOLMatrix::compact() {
#pragma omp parallel for
        for (int i = 0; i < nrows(); ++i)
            rows[i].compact();
    }

    void LOLMatrix::drop_zero_columns() {
        compact();

        // mark non-zero columns
        int ncols = count_columns();
        std::vector<bool> zero_cols(ncols, true);

        for (int i = 0; i < nrows(); ++i) // vector<bool> isn't thread safe
            for (const SparseElement &e : rows[i])
                zero_cols[e.index] = false;

        // remap the columns
        std::vector<int> old2new(ncols);
        int new_col_count = 0;
        for (int col = 0; col < ncols; ++col)
            if (!zero_cols[col])
                old2new[col] = new_col_count++;

        // apply the map
#pragma omp parallel for
        for (int i = 0; i < nrows(); ++i)
            for (SparseElement &e : rows[i])
                e.index = old2new[e.index];
    }

    int LOLMatrix::count_columns() const {
        int max_index = -1;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(max:max_index)
#endif
        for (int i = 0; i < nrows(); ++i)
            for (const SparseElement &e : rows[i])
                max_index = std::max(max_index, e.index);
        return max_index + 1;
    }

    int LOLMatrix::count_nnz() const {
        int cnt = 0;
#pragma omp parallel for reduction(+:cnt)
        for (int i = 0; i < nrows(); ++i)
            cnt += rows[i].size();
        return cnt;
    }

    CRSMatrix LOLMatrix::to_crs() {
        compact();
        int nnz = count_nnz();
        CRSMatrix m = {
            std::vector(nnz, SparseElement(0, 0.)),
            std::vector(nrows()+1, 0)
        };
        for (int i=0; i<nrows(); ++i) {
            m.offset[i+1] = m.offset[i];
            for (const SparseElement &e : rows[i].data)
                m.mat[m.offset[i+1]++] = e;
        }
        return m;
    }
}

