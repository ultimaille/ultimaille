#ifndef __RANGE_H__
#define __RANGE_H__

#include <tuple>
#include <utility>
#include <iterator>

namespace UM {
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
		    T iterable;
	//      typedef decltype(std::begin(std::declval<T>())) iterator_type;
	//      iterator_type iter;
		    bool operator!=(const iterator& rhs) const { return i != rhs.i; }
		    void operator++() { ++i; /*++iter;*/ }
		    auto operator*() const { return std::tie(i, *(iterable.begin()+i)); }
		    auto begin() { return iterator{0, std::forward<T>(iterable)}; }
		    auto end()   { return iterator{std::end(iterable)-std::begin(iterable), std::forward<T>(iterable)}; }
		};
		return iterator{0, std::forward<T>(iterable)};
	}

	#if 0
	template <typename... T> struct zip_helper {
		struct iterator /* : std::iterator<std::forward_iterator_tag, std::tuple<decltype(*std::declval<T>().begin())...>>*/ {
		      using value_ref_tuple_t = std::tuple<typename std::iterator_traits<T>::reference...>;

		        std::tuple<decltype(std::declval<T>().begin())...> iters_;

		        template <std::size_t... I>
		            auto deref(std::index_sequence<I...>) const {
		                return 0;
	//typename iterator::value_type{*std::get<I>(iters_)...};
		            }

		        template <std::size_t... I>
		            void increment(std::index_sequence<I...>) {
		                auto l = {(++std::get<I>(iters_), 0)...};
		            }

		        explicit iterator(decltype(iters_) iters) : iters_{std::move(iters)} {}

		        iterator& operator++() {
		            increment(std::index_sequence_for<T...>{});
		            return *this;
		        }

		        iterator operator++(int) {
		            auto saved{*this};
		            increment(std::index_sequence_for<T...>{});
		            return saved;
		        }

		        bool operator!=(const iterator& other) const {
		            return iters_ != other.iters_;
		        }

		        auto operator*() const { return deref(std::index_sequence_for<T...>{}); }
		};

		zip_helper(T&... seqs)
		    : begin_{std::make_tuple(seqs.begin()...)},
		    end_{std::make_tuple(seqs.end()...)} {}

		iterator begin() const { return begin_; }
		iterator end() const { return end_; }

		iterator begin_;
		iterator end_;
	};

	template <typename... T> constexpr auto zip(T && ... seqs) {
		return zip_helper<T...>{seqs...};
	}

	#endif
}
#endif // __RANGE_H__

