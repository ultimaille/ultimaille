#ifndef __HILBERT_SORT_H__
#define __HILBERT_SORT_H__

#include <vector>
#include "vec.h"
#include "pointset.h"


namespace UM {
    //    _  _ _ _ _             _     ___          _
    //   | || (_) | |__  ___ _ _| |_  / __| ___ _ _| |_
    //   | __ | | | '_ \/ -_) '_|  _| \__ \/ _ \ '_|  _|
    //   |_||_|_|_|_.__/\___|_|  \__| |___/\___/_|  \__|
    //
    // TODO: vec3* -> PointSet
    struct HilbertSort {
        typedef std::vector<int>::iterator Iter;

        HilbertSort(vec3* p_pts, int p_nb_pts) : pts(p_pts), nb_pts(p_nb_pts) {}

        void apply(std::vector<int>& ind, int begin, int end) {
            assert(begin >= 0 && begin <= ind.size());
            assert(end   >= 0 && end   <= ind.size());
            hilbert_sort<0, false, false, false>(ind.begin() + begin, ind.begin() + end);
        }

        // defines 6 order relations on vec3 (3 axis * 2 directions)
        template <int AX, bool DIR>
            struct CmpHilbert {
                CmpHilbert(vec3* p_pts, int p_nb_pts) {
                    pts = p_pts;
                    nb_pts = p_nb_pts;
                }
                bool operator() (int i, int j) {
                    assert(i  >= 0 &&  i < nb_pts);
                    assert(j  >= 0 &&  j < nb_pts);
                    assert(AX >= 0 && AX < 3);
                    return DIR ?
                        (pts[i][AX] < pts[j][AX]) :
                        (pts[i][AX] > pts[j][AX]);
                }
                int nb_pts;
                vec3* pts;
            };

        // find the median w.r.t "Cmp" order relation
        template <class Cmp>
            Iter hilbert_split(Iter begin, Iter end, const Cmp& cmp) {
                ogf_debug_assert(end >= begin);
                if (begin == end) return begin;

                Iter middle = begin + (end - begin) / 2;
                std::nth_element(begin, middle, end, cmp);

                return middle;
            }

        // do the job
        template <int x, bool upx, bool upy, bool upz>
            void hilbert_sort(Iter begin, Iter end) {
                const int y = (x + 1) % 3, z = (x + 2) % 3;
                if (end - begin <= 2) return;

                Iter m0 = begin, m8 = end;
                Iter m4 = hilbert_split(m0, m8, CmpHilbert<x,  upx>(pts, nb_pts));
                Iter m2 = hilbert_split(m0, m4, CmpHilbert<y,  upy>(pts, nb_pts));
                Iter m1 = hilbert_split(m0, m2, CmpHilbert<z,  upz>(pts, nb_pts));
                Iter m3 = hilbert_split(m2, m4, CmpHilbert<z, !upz>(pts, nb_pts));
                Iter m6 = hilbert_split(m4, m8, CmpHilbert<y, !upy>(pts, nb_pts));
                Iter m5 = hilbert_split(m4, m6, CmpHilbert<z,  upz>(pts, nb_pts));
                Iter m7 = hilbert_split(m6, m8, CmpHilbert<z, !upz>(pts, nb_pts));

                hilbert_sort<z,  upz,  upx,  upy>(m0, m1);
                hilbert_sort<y,  upy,  upz,  upx>(m1, m2);
                hilbert_sort<y,  upy,  upz,  upx>(m2, m3);
                hilbert_sort<x,  upx, !upy, !upz>(m3, m4);
                hilbert_sort<x,  upx, !upy, !upz>(m4, m5);
                hilbert_sort<y, !upy,  upz, !upx>(m5, m6);
                hilbert_sort<y, !upy,  upz, !upx>(m6, m7);
                hilbert_sort<z, !upz, !upx,  upy>(m7, m8);
            }

        vec3* pts;
        int nb_pts;
    };
}

#endif //__HILBERT_SORT_H__

