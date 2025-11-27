#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <iostream>
#include <cmath>
#include <vector>

namespace UM {
    struct SparseElement {
        SparseElement() = default;
        SparseElement(const SparseElement  &) = default;
        SparseElement(      SparseElement &&) = default;
        SparseElement& operator=(const SparseElement  &) = default;
        SparseElement& operator=(      SparseElement &&) = default;
        SparseElement(int index, double value) : index(index), value(value) {}

        inline bool is_null() const { return std::abs(value) < TOL; }

        int index = 0;
        double value = 0.;
        static constexpr double TOL = 1e-10; // tolerance for the element to be considered null
    };

    inline SparseElement operator*(const SparseElement& t, double c) {
        return { t.index, c * t.value };
    }

    inline SparseElement operator*(double c, const SparseElement& t) {
        return t*c;
    }

    inline SparseElement operator-(const SparseElement& t) {
        return { t.index, -t.value };
    }
    std::ostream& operator<<(std::ostream& out, const SparseElement& e);

    //////////////////////////////////////////////////////////////////////////////////////////////////

    struct SparseVector {
        SparseVector(SparseElement e) : data{e} {}
        SparseVector(std::initializer_list<SparseElement> l) : data{l} { compact(); }
        SparseVector(std::vector<SparseElement> &&v) : data{std::move(v)} { compact(); }

        // sort by index, sum same-index terms, remove near-zero entries
        void compact();

        inline SparseVector& operator+=(const SparseVector& other) {
            data.insert(data.end(), other.begin(), other.end());
            compact();
            return *this;
        }

        inline SparseVector& operator-=(const SparseVector& other) {
            data.reserve(size() + other.size());
            for (const SparseElement &e : other)
                data.push_back(-e);
            compact();
            return *this;
        }

        inline SparseVector& operator*=(double x) {
            for (SparseElement &e : data)
                e.value *= x;
            compact();
            return *this;
        }

        // basic wrappers around the container
        SparseVector() = default;
        SparseVector(const SparseVector  &) = default;
        SparseVector(      SparseVector &&) = default;
        SparseVector& operator=(const SparseVector  &) = default;
        SparseVector& operator=(      SparseVector &&) = default;

        inline int size() const { return static_cast<int>(data.size()); }
        inline bool empty() const { return data.empty(); }
        inline       SparseElement& front()       { return data.front(); }
        inline const SparseElement& front() const { return data.front(); }
        inline       SparseElement& back()       { return data.back(); }
        inline const SparseElement& back() const { return data.back(); }
        inline       SparseElement& operator[](int i)       { return data[i]; }
        inline const SparseElement& operator[](int i) const { return data[i]; }

        inline std::vector<SparseElement>::iterator begin() { return data.begin(); }
        inline std::vector<SparseElement>::iterator end()   { return data.end();   }
        inline std::vector<SparseElement>::const_iterator begin() const { return data.begin(); }
        inline std::vector<SparseElement>::const_iterator end()   const { return data.end();   }

        std::vector<SparseElement> data = {};
    };

    std::ostream& operator<<(std::ostream& out, const SparseVector& v);

    inline SparseVector operator+(const SparseVector& a, const SparseVector& b) {
        SparseVector c;
        c.data.reserve(a.size() + b.size());
        c.data.insert(c.end(), a.begin(), a.end());
        c.data.insert(c.end(), b.begin(), b.end());
        c.compact();
        return c;
    }

    inline SparseVector operator-(const SparseVector& a, const SparseVector& b) {
        SparseVector c;
        c.data.reserve(a.size() + b.size());
        c.data.insert(c.end(), a.begin(), a.end());
        for (const SparseElement &e : b)
            c.data.push_back(-e);
        c.compact();
        return c;
    }

    inline SparseVector operator*(const SparseVector& v, double a) {
        SparseVector u = v;
        u *= a;
        return u;
    }

    inline SparseVector operator*(double a, const SparseVector& v) {
        return v * a;
    }

    inline double operator*(const SparseVector& v1, std::vector<double>& v2) {
        double dot = 0;
        for (const SparseElement &e : v1)
            dot += e.value * v2[e.index];
        return dot;
    }

}

#endif //__VECTOR_H__

