#ifndef __HILBERT_SORT_H__
#define __HILBERT_SORT_H__

#include <vector>
#include <cassert>
#include <algorithm>
#include "ultimaille/algebra/vec.h"

namespace UM {
    //    _  _ _ _ _             _     ___          _
    //   | || (_) | |__  ___ _ _| |_  / __| ___ _ _| |_
    //   | __ | | | '_ \/ -_) '_|  _| \__ \/ _ \ '_|  _|
    //   |_||_|_|_|_.__/\___|_|  \__| |___/\___/_|  \__|
    //
    struct HilbertSort {
        HilbertSort(const std::vector<vec3> &data) : pts(data) {}

        void apply(std::vector<int>& ind) const {
            apply(ind, 0, pts.size());
        }

        void apply(std::vector<int>& ind, int begin, int end) const {
            assert(begin >= 0 && begin <= (int)ind.size());
            assert(end   >= 0 && end   <= (int)ind.size());
            hilbert_sort<0, false, false, false>(ind.begin() + begin, ind.begin() + end);
        }

        // defines 6 order relations on vec3 (3 axis * 2 directions)
        template <int AX, bool DIR>
            struct CmpHilbert {
                CmpHilbert(const std::vector<vec3> &data) : pts(data) {}
                bool operator() (int i, int j) const {
                    assert(i  >= 0 &&  i < (int)pts.size());
                    assert(j  >= 0 &&  j < (int)pts.size());
                    assert(AX >= 0 && AX < 3);
                    return DIR ?
                        (pts[i][AX] < pts[j][AX]) :
                        (pts[i][AX] > pts[j][AX]);
                }
                const std::vector<vec3> &pts;
            };

        // find the median w.r.t "Cmp" order relation
        template <class Cmp>
            std::vector<int>::iterator hilbert_split(std::vector<int>::iterator begin, std::vector<int>::iterator end, const Cmp& cmp) const {
                assert(end >= begin);
                if (begin == end) return begin;
                std::vector<int>::iterator middle = begin + (end - begin) / 2;
                std::nth_element(begin, middle, end, cmp);
                return middle;
            }

        // do the job
        template <int x, bool upx, bool upy, bool upz>
            void hilbert_sort(std::vector<int>::iterator begin, std::vector<int>::iterator end) const {
                constexpr int y = (x + 1)%3, z = (x + 2)%3;
                if (end - begin <= 2) return;

                std::vector<int>::iterator m0 = begin, m8 = end;
                std::vector<int>::iterator m4 = hilbert_split(m0, m8, CmpHilbert<x,  upx>(pts));
                std::vector<int>::iterator m2 = hilbert_split(m0, m4, CmpHilbert<y,  upy>(pts));
                std::vector<int>::iterator m1 = hilbert_split(m0, m2, CmpHilbert<z,  upz>(pts));
                std::vector<int>::iterator m3 = hilbert_split(m2, m4, CmpHilbert<z, !upz>(pts));
                std::vector<int>::iterator m6 = hilbert_split(m4, m8, CmpHilbert<y, !upy>(pts));
                std::vector<int>::iterator m5 = hilbert_split(m4, m6, CmpHilbert<z,  upz>(pts));
                std::vector<int>::iterator m7 = hilbert_split(m6, m8, CmpHilbert<z, !upz>(pts));

                hilbert_sort<z,  upz,  upx,  upy>(m0, m1);
                hilbert_sort<y,  upy,  upz,  upx>(m1, m2);
                hilbert_sort<y,  upy,  upz,  upx>(m2, m3);
                hilbert_sort<x,  upx, !upy, !upz>(m3, m4);
                hilbert_sort<x,  upx, !upy, !upz>(m4, m5);
                hilbert_sort<y, !upy,  upz, !upx>(m5, m6);
                hilbert_sort<y, !upy,  upz, !upx>(m6, m7);
                hilbert_sort<z, !upz, !upx,  upy>(m7, m8);
            }

        const std::vector<vec3> &pts;
    };
}

#endif //__HILBERT_SORT_H__

