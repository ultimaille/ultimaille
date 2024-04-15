#ifndef __LINEXPR_H__
#define __LINEXPR_H__

#include <initializer_list>
#include <iostream>
#include <vector>

namespace UM {
    namespace Linear {

        // Sparse representation of a linear expression b + \sigma_i a_i x_i.
        // A linear expression consists of a constant term, plus a list of coefficient-variable pairs that capture the linear terms.
        // Linear expressions are used to build linear objective and constraints.
        // They are temporary objects that typically have short lifespans.

        struct LinExpr {
            struct Term {
                Term() {}
                Term(int i, double a) : i(i), a(a) {}
                Term(double a) : a(a) {}

                inline bool is_null() const { return std::abs(a) < TOL; }

                int i = -1; // -1 means constant term (b)
                double a = 0.;
                static constexpr double TOL = 1e-10;
            };

            LinExpr() = default;
            LinExpr(Term c) : data{c} {}
            LinExpr(double a) : data{{a}} {}
            LinExpr(std::initializer_list<Term> c) : data{c} { compact(); }

            // sort (by i) and aggreate terms, remove near-zero entries
            void compact();

            // N.B. expressions must be compacted before doing any arithmetics
            LinExpr& operator+=(const LinExpr& e);
            LinExpr& operator-=(const LinExpr& e);
            LinExpr& operator*=(double a);

            // basic wrappers around the container
            inline int size() const { return data.size(); }
            inline bool empty() const { return data.empty(); }
            inline       Term& back()       { return data.back(); }
            inline const Term& back() const { return data.back(); }
            inline       Term& operator[](int i)       { return data[i]; }
            inline const Term& operator[](int i) const { return data[i]; }
            std::vector<Term>::iterator begin() { return data.begin(); }
            std::vector<Term>::iterator end()   { return data.end();   }
            std::vector<Term>::const_iterator begin() const { return data.begin(); }
            std::vector<Term>::const_iterator end()   const { return data.end();   }

            std::vector<Term> data = {}; // TODO initial capacity?
        };

        inline LinExpr X(int i) { return {{i, 1.}}; }

        std::ostream& operator<<(std::ostream& os, const LinExpr& le);

        inline bool operator<(const LinExpr::Term& t1, const LinExpr::Term& t2) {
            return t1.i < t2.i;
        }

        inline bool operator!=(const LinExpr::Term& t1, const LinExpr::Term& t2) {
            return t1.i != t2.i || std::abs(t1.a - t2.a) > LinExpr::Term::TOL;
        }

        inline bool operator==(const LinExpr::Term& t1, const LinExpr::Term& t2) {
            return !(t1 != t2);
        }

        inline LinExpr::Term operator-(const LinExpr::Term &t) {
            return { t.i, -t.a };
        }

        inline LinExpr::Term operator*(const LinExpr::Term& t, double c) {
            return { t.i, c*t.a };
        }

        inline LinExpr::Term operator*(double c, const LinExpr::Term& t) {
            return t*c;
        }

        LinExpr operator+(const LinExpr& a, const LinExpr& b);
        LinExpr operator-(const LinExpr& a, const LinExpr& b);
        LinExpr operator*(const LinExpr& le, double c);
        inline LinExpr operator*(double c, const LinExpr& le) {
            return le * c;
        }

        inline LinExpr& LinExpr::operator+=(const LinExpr& e) { // TODO measure whether we need to define += or + first
            return *this = LinExpr(*this) + e;
        }

        inline LinExpr& LinExpr::operator-=(const LinExpr& e) {
            return *this = LinExpr(*this) - e;
        }

        inline LinExpr& LinExpr::operator*=(double c) {
            return *this = LinExpr(*this) * c;
        }
    }
}

#endif //__LINEXPR_H__

