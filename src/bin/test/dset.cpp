#undef NDEBUG
#include <cassert>
#include <iostream>
#include <vector>
#include <utility>
#include <ctime>

#include "ultimaille/disjointset.h"

int main() {
    const int n=1017, maxnsets=13;
    std::srand((unsigned int)std::time(nullptr));
    std::vector<int> values(n);
    for (auto &v : values) v = (std::rand()%maxnsets);

    // create the list of pairs
    std::vector<std::pair<int,int> > perm;
    for (int i=0; i<n; ++i)
        for (int j=i+1; j<n; ++j)
            if (values[i]==values[j])
                perm.push_back(std::make_pair(i,j));

    // randomize this list
    for (int j=perm.size(); j>0; --j)
        std::swap(perm[j-1], perm[std::rand()%j]);

    {
        DisjointSet dSet(n);
        for (auto const &v : perm)
            dSet.merge(v.first, v.second);

        std::vector<int> indices;
        int nsets = dSet.get_sets_id(indices);

        assert(indices.size()==n);
        assert(nsets<=maxnsets);

        for (int i=0; i<n; ++i)
            for (int j=i+1; j<n; ++j)
                assert((values[i]==values[j])==(indices[i]==indices[j]));
    }

    {
       for (auto &v : values) v *= (std::rand()%2 ? 1 : -1);
       DisjointSetWithSign dSet(n);
       for (auto const &v : perm)
           dSet.merge(v.first, v.second, values[v.first]==values[v.second]);


/*
   GEO::vector<size_t> indices;
   dSet.getSetsId(indices);
   if (indices.size()!=values.size())
    return false;
   auto const &signs=dSet.getSameSigns();
   if (indices.size()!=signs.size())
    return false;
   for (size_t i=0; i<values.size(); ++i) {
    for (size_t j=i+1; j<values.size(); ++j) {
     if ((std::abs(values[i])==std::abs(values[j])) != (indices[index_t(i)]==indices[index_t(j)]))
      return false;
     if (std::abs(values[i])!=std::abs(values[j]))
      continue;
     if ((values[i]==values[j])!=(signs[i]==signs[j]))
      return false;
    }
   }
*/

    }

    return 0;
}

