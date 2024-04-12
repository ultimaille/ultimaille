#include "linexpr.h"
#include <algorithm>

namespace UM {
    namespace Linear {
        void LinExpr::compact() {
            std::sort(data.begin(), data.end());
            int s = 0;
            for (int i = 0; i < (int)data.size(); ) {
                data[s] = data[i++];
                while (i < (int)data.size() && data[i].i == data[s].i)
                    data[s].a += data[i++].a;

                if (std::abs(data[s].a) > LinExpr::Term::TOL)
                    s++;
            }
            data.resize(s);
        }

        LinExpr operator+(const LinExpr& le1, const LinExpr& le2) {
            LinExpr result;
            int i = 0, j = 0;
            while (i < (int)le1.data.size() && j < (int)le2.data.size()) {
                if (le1.data[i].i == le2.data[j].i) {
                    const double a = le1.data[i].a + le2.data[j].a;
                    if (std::abs(a) > LinExpr::Term::TOL)
                        result.data.emplace_back(le1.data[i].i, a);
                    ++i; ++j;
                } else if (le1.data[i].i < le2.data[j].i) {
                    result.data.push_back(le1.data[i++]);
                } else {
                    result.data.push_back(le2.data[j++]);
                }
            }
            result.data.reserve(result.data.size() + le1.data.size() + le2.data.size() - (i + j));
            result.data.insert(result.data.end(), le1.data.begin() + i, le1.data.end());
            result.data.insert(result.data.end(), le2.data.begin() + j, le2.data.end());
            return result;
        }

        LinExpr operator-(const LinExpr& le1, const LinExpr& le2) {
            LinExpr result;
            int i = 0, j = 0;
            while (i < (int)le1.data.size() && j < (int)le2.data.size()) {
                if (le1.data[i].i == le2.data[j].i) {
                    const double a = le1.data[i].a + le2.data[j].a;
                    if (std::abs(a) > LinExpr::Term::TOL)
                        result.data.emplace_back(le1.data[i].i, a);
                    ++i; ++j;
                } else if (le1.data[i].i < le2.data[j].i) {
                    result.data.push_back(le1.data[i++]);
                } else {
                    result.data.push_back(-le2.data[j++]);
                }
            }
            result.data.reserve(result.data.size() + le1.data.size() + le2.data.size() - (i + j));
            result.data.insert(result.data.end(), le1.data.begin() + i, le1.data.end());
            while (j < (int)le2.data.size()) result.data.push_back(-le2.data[j++]);
            return result;
        }

        LinExpr operator*(const LinExpr& le, double c) {
            if (c == 0.) return {};
            LinExpr result;
            result.data.reserve(le.data.size());
            double eps = LinExpr::Term::TOL / std::abs(c);
            for (const LinExpr::Term &t : le.data)
                if (std::abs(t.a) > eps)
                    result.data.push_back(c * t);
            return result;
        }

        std::ostream& operator<<(std::ostream& os, const LinExpr& le) {
            for (auto &t : le.data) {
                if (t.a < 0)
                    os << " - ";
                else
                    os << " + ";

                if (t.i < 0)
                    os << std::abs(t.a);
                else {
                    if (std::abs(t.a) != 1)
                        os << std::abs(t.a) << "*";
                    os << "X_" << t.i;
                }
            }
            return os;
        }
    }
}

