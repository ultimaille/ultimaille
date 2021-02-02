#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#undef NDEBUG
#include <cassert>

#include <ultimaille/permutation.h>

using namespace UM;

int main() {
    int n = 1277;
    std::vector<int> destinations(n);
    std::iota(destinations.begin(), destinations.end(), 0);
    std::vector<int> data = destinations;

    std::shuffle(destinations.begin(), destinations.end(), std::mt19937{std::random_device{}()});

    Permutation p(destinations);
    p.apply(data);

    for (int i=0; i<n; i++)
        assert(destinations[i]==data[i]);

    return 0;
}

