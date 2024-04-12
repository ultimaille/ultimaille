#ifndef __RBUILDER_H__
#define __RBUILDER_H__

#include "linexpr.h"
#include "sparsematrix.h"
#include <vector>

namespace UM {
    namespace Linear {
        struct ReductionBuilder {

            CRSMatrix get_matrix();

            /*
            private:
                struct ConstraintMatrix { // nvars x (nvars + 1) matrix (+1 comes for the affine constant)
                    ConstraintMatrix(int nvars) : nvars(nvars), constraints(nvars) { // identity matrix
                        for (int i=0; i<nvars; i++)
                            constraints[i] = X(i);
                    }

                    int nvars;
                    std::vector<LinExpr> constraints;
                };

                ConstraintMatrix M;
                */
        };

#if 0
        struct ConstraintMatrix {
		void drop_zero_columns() {
			for (LinExpr &le : constraints)
				le.compact();

			int nvars() = count_variables(); // does not count the affine constant

			// mark non-zero columns
			std::vector<bool> zero_cols(ncols, true);

			zero_cols[0] = false; // do not remove the affine column even if empty
			for (int i = 0; i < (int)constraints.size(); ++i)
				for (const LinExpr::Term &c : constraints[i].data)
					zero_cols[c.i + 1] = false;

			std::vector<int> old2new(ncols);
			int new_col_count = 0;
			for (int col = 0; col < ncols; ++col)
				if (!zero_cols[col])
					old2new[col] = new_col_count++;

			// compact
			for (int row = 0; row < (int)constraints.size(); ++row) {
				for (const LinExpr::Term &c : constraints[i].data) {
					if 
				}

				int writeIndex = 0;
				for (size_t i = 0; i < columns[row].size(); ++i) {
					int col = columns[row][i];
					if (!zero_cols[col]) {
						columns[row][writeIndex] = columnRemap[col];
						values[row][writeIndex] = values[row][i];
						++writeIndex;
					}
				}
				columns[row].resize(writeIndex);
				values[row].resize(writeIndex);
			}
		}

		int count_variables() const { // N.B. does not account for the constant!
			int max_index = -1;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(max:max_index)
#endif
			for (int i = 0; i < (int)constraints.size(); ++i)
				for (const LinExpr::Term &c : constraints[i].data)
					max_index = std::max(max_index, c.i);
			return max_index + 1;
		}

            std::vector<LinExpr> constraints = {};
        };

//      inline ConstraintMatrix cm_identity(int n) {
//          ConstraintMatrix m(n, n);
//          for (int i=0; i<n; i++)
//              m.rows[i] = X(i);
//      }
#endif
    }
}

#endif //__RBUILDER_H__

