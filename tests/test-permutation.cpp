#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("Permutation", "[permutation]") {
    int n = 1775;
    Permutation perm(n);
    CHECK( perm.is_valid() );

    std::shuffle(perm.begin(), perm.end(), std::mt19937{std::random_device{}()});
    CHECK( perm.is_valid() );

    std::vector<int> data(n);
    std::iota(data.begin(), data.end(), 0);

    perm.apply(data);
    for (int i=0; i<n; i++)
        CHECK( perm[i]==data[i] );

    perm.apply_reverse_lin(data);
    for (int i=0; i<n; i++)
        CHECK( data[i]==i );

    perm.apply_lin(data);
    for (int i=0; i<n; i++)
        CHECK( perm[i]==data[i] );

    perm.apply_reverse(data);
    for (int i=0; i<n; i++)
        CHECK( data[i]==i );


    {
        std::vector<char> A = {'a', 'b', 'c', 'd', 'e'};
        std::vector<char> B = {'e', 'd', 'c', 'a', 'b'};
        Permutation perm(5);
        perm.ind = {4, 3, 2, 0, 1};
        perm.apply(A);
        for (int i=0; i<5; i++)
            CHECK( A[i]==B[i] );
//      for (char c : A)
//          std::cerr << c << " ";
//      std::cerr << std::endl;
    }

    {
        std::vector<char> A = {'a', 'b', 'c', 'd', 'e'};
        std::vector<char> B = {'e', 'd', 'c', 'a', 'b'};
        Permutation perm(5);
        perm.ind = {4, 3, 2, 0, 1};
        perm.apply_reverse(B);
        for (int i=0; i<5; i++)
            CHECK( A[i]==B[i] );
//      for (char c : B)
//          std::cerr << c << " ";
//      std::cerr << std::endl;
    }
}

