#include "covariance.h"

namespace UM {
    PointSetCovariance operator+(PointSetCovariance& A, PointSetCovariance& B) {
        PointSetCovariance res;
        res.n = A.n + B.n;
        res.center = (A.center * A.n + B.center * B.n) /res.n;
        vec3 delta = B.center - A.center;
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                res.cov[i][j] = (A.cov[i][j]*A.n + B.cov[i][j]*B.n + (delta[i] * delta[j] * A.n*B.n)/res.n)/res.n;
        return res;
    }
}

