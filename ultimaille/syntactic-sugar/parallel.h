#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#include <atomic>

namespace UM {

#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
    initializer(omp_priv = std::vector<double>(omp_orig.size(), 0))
#endif

    struct SpinLock {
        SpinLock() : flag{ATOMIC_FLAG_INIT} {}

        void lock() {
            while (flag.test_and_set(std::memory_order_acquire));
        }

        void unlock() {
            flag.clear();
        }
        private:
        std::atomic_flag flag;
    };

}

#endif // __PARALLEL_H__

