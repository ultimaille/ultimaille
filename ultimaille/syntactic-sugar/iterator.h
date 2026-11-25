#ifndef __ITERATOR_H__
#define __ITERATOR_H__

template <typename P, typename M>
auto custom_iterator(M &m, int from, int to) {
    struct iterator {
        P data;
        P& operator*()    { return data; }
        void operator++() { ++data; }
        bool operator!=(iterator& rhs) { return data != rhs.data; }
    };
    struct wrapper {
        M& m;
        int from, to;
        iterator begin() { return {{ m, from }}; }
        iterator end()   { return {{ m, to   }}; }
    };
    return wrapper{ m, from, to };
}

#endif // __ITERATOR_H__

