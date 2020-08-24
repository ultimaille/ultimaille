#ifndef __RANGE_H__
#define __RANGE_H__

#include <tuple>
#include <utility>
#include <iterator>

constexpr auto range(int n) {
    struct iterator {
        int i;
        void operator++() { ++i; }
        bool operator!=(const iterator& rhs) const { return i != rhs.i; }
        const int &operator*() const { return i; }
    };
    struct wrapper {
        int n;
        auto begin() { return iterator{0}; }
        auto end()   { return iterator{n}; }
    };
    return wrapper{n};
}

template <typename T> constexpr auto enumerate(T && iterable) {
    struct iterator {
        int i;
        typedef decltype(std::begin(std::declval<T>())) iterator_type;
        iterator_type iter;
        bool operator!=(const iterator& rhs) const { return iter != rhs.iter; }
        void operator++() { ++i; ++iter; }
        auto operator*() const { return std::tie(i, *iter); }
    };
    struct wrapper {
        T iterable;
        auto begin() { return iterator{0, std::begin(iterable)}; }
        auto end()   { return iterator{0, std::end  (iterable)}; }
    };
    return wrapper{std::forward<T>(iterable)};
}

#endif // __RANGE_H__

