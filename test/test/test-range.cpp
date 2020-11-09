#include <iostream>
#include <cstring>
#include <vector>

#include <ultimaille/all.h>

using namespace UM;

int main() {
    std::cerr << "testing for (int i : range(7)): " << std::endl;
    for (int i : range(7))
        std::cerr << i  << " ";
    std::cerr << std::endl << std::endl;

    std::cerr << "testing for (auto [i,v] : enumerate(array)): " << std::endl;
    std::vector<std::string> array{"a", "b", "c", "d", "e", "f", "g"};
    for (auto [i,v] : enumerate(array))
        std::cerr << "(" << i << ",'" << v << "') ";
    std::cerr << std::endl << std::endl;

    std::cerr << "testing for (auto [u,v] : zip(array_a, array_b)): " << std::endl;
    std::vector<int> array2{3, 1, 4, 1, 5, 9, 2};
    for (auto [u,v] : zip(array2, array))
        std::cerr << "(" << u << ",'" << v << "') ";
    std::cerr << std::endl;

    return 0;
}

