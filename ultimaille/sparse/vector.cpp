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

    std::ostream& operator<<(std::ostream& out, const SparseElement& e) {
        out << "{" << e.index << ", " << e.value << "}";
        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SparseVector& v) {
        out << "{";
        for (int i=0; i<v.size(); i++) {
            out << v[i];
            if (i+1<v.size()) out << ", ";
        }
        out << "}";
        return out;
    }
}

