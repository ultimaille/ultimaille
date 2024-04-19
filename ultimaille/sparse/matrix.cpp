#include "matrix.h"

namespace UM {

    double CRSMatrix::dot(const std::vector<double> &x) const {
        double result = 0;
#pragma omp parallel for reduction(+:result)
        for (int row = 0; row < nrows(); ++row) {
            for (const SparseElement &e : iter_row(row))
                result += e.value * x[e.index];
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
        std::vector<SparseElement> newmat(nnz(), SparseElement{});
        std::vector<int> newoffset(ncols+1, 0);

        // count nnz per column
        for (int row = 0; row < nrows(); row++)
            for (const SparseElement &e : iter_row(row))
                newoffset[1+e.index]++;

        // sum the counts to obtain the column offsets
        for (int col = 0; col < ncols; col++)
            newoffset[col+1] += newoffset[col];

        // fill in the transposed data
        std::vector<int> cur_nnz_per_col(ncols, 0);
        for (int row = 0; row < nrows(); row++)
            for (const SparseElement &e : iter_row(row))
                newmat[ newoffset[e.index] + cur_nnz_per_col[e.index]++ ] = { row, e.value };

        return { newmat, newoffset };
    }

    // given a n-component row vector X, multiply it by a m x n matrix M : Y = X M
    // usually it can be computed as Y = sum_{j=0}^{m-1} (0 0 0 0 0 ... sum_{i=0}^{n-1} Xi Mij ... 0 0 0 0)
    // but we can refactor it as Y = sum_{j=0}^{m-1} sum_{i=0}^{n-1} Xi (0 0 0 0 0 ... Mij ... 0 0 0 0)
    // and further Y = sum_{i=0}^{n-1} Xi sum_{j=0}^{m-1} (0 0 0 0 0 ... Mij ... 0 0 0 0)

    SparseVector operator*(const SparseVector& X, const CRSMatrix& M) {
        std::vector<SparseElement> data;
        for (const SparseElement &Xi : X)
            for (const SparseElement &Mij : M.iter_row(Xi.index))
                data.push_back(Mij * Xi.value);
        return SparseVector(std::move(data));
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
            std::vector<SparseElement>(nnz, SparseElement{}),
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

