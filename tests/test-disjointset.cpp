#include <catch2/catch.hpp>

#include <iostream>
#include <vector>
#include <utility>
#include <ctime>

#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("Test union find", "[DisjointSet]") {
    const int n=10327, maxnsets=394;
    std::srand((unsigned int)std::time(nullptr));
    std::vector<int> values(n);
    for (auto &v : values) v = 1 + std::rand()%maxnsets;

    // create the list of pairs
    std::vector<std::pair<int,int> > perm;
    for (int i=0; i<n; ++i)
        for (int j=i+1; j<n; ++j)
            if (values[i]==values[j])
                perm.push_back(std::make_pair(i,j));

    // randomize this list
    for (int j=perm.size(); j>0; --j)
        std::swap(perm[j-1], perm[std::rand()%j]);

    { // test the regular disjoint set
        DisjointSet dSet(n);
        for (auto const &v : perm)
            dSet.merge(v.first, v.second);

        std::vector<int> indices;
        int nsets = dSet.get_sets_id(indices);

        CHECK( indices.size()==n );
        CHECK( nsets<=maxnsets );

        for (int i=0; i<n; ++i)
            for (int j=i+1; j<n; ++j)
                CHECK( (values[i]==values[j]) == (indices[i]==indices[j]) );
    }
}

