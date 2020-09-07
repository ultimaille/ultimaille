#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"
#include "ultimaille/range.h"
#include <cstdio>

#include <algorithm>

#include <list>

#include "ZipIterator.hpp" //Header only ;)

/*
inline auto facets(Surface &m) {
    return range(m.nfacets());
}

inline auto corners(Surface &m) {
    return range(m.ncorners());
}

struct facet_vertices {
    Surface &m_;
    const int facet_;
    facet_vertices(Surface &m, const int facet) : m_(m), facet_(facet) {}

    struct iterator {
        void operator++() { ++value_; }
        int &operator*() const { return m_.vert(facet_, value_); }
        bool operator!=(const iterator& rhs) const { return value_ != rhs.value_; }
        Surface &m_;
        const int facet_;
        int value_;
    };

    iterator begin() const { return iterator{m_, facet_, 0};    }
    iterator end()   const { return iterator{m_, facet_, m_.facet_size(facet_)}; }
};
*/

#if 0


#include <utility>

//namespace impl {

    template <typename Iter, typename... Iters>
        class Zip_iterator {
            public:
                using value_type = std::tuple<const typename Iter::value_type&, const typename Iters::value_type&...>;

                Zip_iterator(const Iter &head, const Iters&... tail)
                    : head_(head), tail_(tail...) { }

                value_type operator*() const {
                    return std::tuple_cat(std::tuple<const typename Iter::value_type&>(*head_), *tail_);
                }

                Zip_iterator& operator++() {
                    ++head_; ++tail_;
                    return *this;
                }

                bool operator==(const Zip_iterator &rhs) const {
                    return head_ == rhs.head_ && tail_ == rhs.tail_;
                }

                bool operator!=(const Zip_iterator &rhs) const {
                    return !(*this == rhs);
                }

            private:
                Iter head_;
                Zip_iterator<Iters...> tail_;
        };

    template <typename Iter>
        class Zip_iterator<Iter> {
            public:
                using value_type = std::tuple<const typename Iter::value_type&>;

                Zip_iterator(const Iter &head) : head_(head) { }

                value_type operator*() const {
                    return value_type(*head_);
                }

                Zip_iterator& operator++() { ++head_; return *this; }

                bool operator==(const Zip_iterator &rhs) const { return head_ == rhs.head_; }

                bool operator!=(const Zip_iterator &rhs) const { return !(*this == rhs); }

            private:
                Iter head_;
        };

//}  // namespace impl

template <typename Iter>
class seq {
public:
  using iterator = Iter;
  seq(const Iter &begin, const Iter &end) : begin_(begin), end_(end) { }
  iterator begin() const { return begin_; }
  iterator end() const { return end_; }
private:
  Iter begin_, end_;
};

/* WARNING: Undefined behavior if iterator lengths are different.
 */
template <typename... Seqs>
seq<Zip_iterator<typename Seqs::iterator...>>
Zip(const Seqs&... seqs) {
  return seq<Zip_iterator<typename Seqs::iterator...>>(
      Zip_iterator<typename Seqs::iterator...>(std::begin(seqs)...),
      Zip_iterator<typename Seqs::iterator...>(std::end(seqs)...));
}

#endif

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }

    Polygons pm;
    SurfaceAttributes attributes = read_geogram(argv[1], pm);

/*

    Triangles tri;
    pm.extract_triangles(tri);
    assert(tri.nfacets() == pm.nfacets());

    Quads quads;
    pm.extract_quads(quads);
//  assert(quads.nfacets() == pm.nfacets());
//  std::cerr << quads.points.use_count() << std::endl;

//  PointAttribute<int> prand("rand", attributes, tri);
//  FacetAttribute<int> fid("id", attributes, tri);
//  CornerAttribute<int> cid("id", attributes, tri);

    FacetAttribute<int> nonex("nonexisting", attributes, tri);
    for (int i=0; i<pm.nfacets(); i++)
        nonex[i] = rand()%1980;


*/
//  for (int f : facets(pm))
//      for (int &v : facet_vertices(pm, f))
//          v = rand()%999;

/*
    for (int f : facets(pm))
        for (auto [lv, v] : enumerate(facet_vertices(pm, f))) {
            std::cerr << lv << " " << v << std::endl;
        }
*/
    write_geogram("read_test.geogram", pm, attributes);

    std::vector<int> one{{111, 1, 11}};
    const std::list<float> three{{333, 3, 33}};
    std::vector<int> four{{444, 4, 44, -11111}};
    for (auto [i,v] : enumerate(one)) {
        std::cerr << i << " " << v << std::endl;
    }
    return 0;
    /*
    for (auto a : Zip(one, three, four)) {
        std::printf("%d %f %d\n", std::get<0>(a),
                std::get<const float&>(a), std::get<2>(a));
    }
    */

    for (auto [v1, v2, v3]: Zip(one, three, four)) std::cout << v1++ << " " << v2  << " " << v3++ << std::endl;
    std::cout << std::endl;

//  auto zip = Zip(one, three, four);
//    std::sort(zip.begin(), zip.end());

    for (auto [v1, v2, v3]: Zip(one, three, four)) std::cout << v1 << " " << v2 << " " << v3 << std::endl;
    std::cout << std::endl;


    return 0;
}

