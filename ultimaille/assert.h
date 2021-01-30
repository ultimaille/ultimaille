#ifndef __ASSERT_H__
#define __ASSERT_H__

#include <iostream>
#include <cstdlib>

namespace UM {
    struct source_location {
        const char* file_name;
        unsigned line_number;
        const char* function_name;
    };

#if defined(__linux__)
#define CUR_SOURCE_LOCATION source_location{__FILE__, __LINE__, __PRETTY_FUNCTION__}
#else
#define CUR_SOURCE_LOCATION source_location{__FILE__, __LINE__, __FUNCSIG__}
#endif

    void release_assert(bool expr, const source_location& loc, const char* expression) {
        if (expr) return;
        std::cerr << loc.file_name << ":" << loc.line_number << ": " << loc.function_name << ": Assertion `" << expression << "` failed." << std::endl;
        std::abort();
    }

#define um_assert(Expr) release_assert(Expr, CUR_SOURCE_LOCATION, #Expr)
}

#endif //__ASSERT_H__

