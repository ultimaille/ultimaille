#include <iostream>
#include <cstdlib>

#include "ultimaille/syntactic-sugar/assert.h"

namespace UM {
    void release_assert(bool expr, const source_location& loc, const char* expression) {
        if (expr) return;
        std::cerr << loc.file_name << ":" << loc.line_number << ": " << loc.function_name << ": Assertion `" << expression << "` failed." << std::endl;
        std::abort();
    }
}

