#include "vector.h"
#include <algorithm>

namespace UM {
    void SparseVector::compact() {
        std::sort(begin(), end(),  [](const SparseElement& a, const SparseElement& b) { return a.index < b.index; });
        int s = 0;
        for (int i = 0; i < size(); ) {
            data[s] = data[i++];
            while (i < size() && data[i].index == data[s].index)
                data[s].value += data[i++].value;

            if (!data[s].is_null())
                s++;
        }
        data.resize(s);
    }
}

