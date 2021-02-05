#ifndef __ASSERT_H__
#define __ASSERT_H__

namespace UM {
    struct source_location {
        const char* file_name;
        unsigned line_number;
        const char* function_name;
    };

    void release_assert(bool expr, const source_location& loc, const char* expression);

#if defined(__GNUC__)
#define CUR_SOURCE_LOCATION source_location{__FILE__, __LINE__, __PRETTY_FUNCTION__}
#else
#define CUR_SOURCE_LOCATION source_location{__FILE__, __LINE__, __FUNCSIG__}
#endif

#define um_assert(Expr) release_assert(Expr, CUR_SOURCE_LOCATION, #Expr)
}

#endif //__ASSERT_H__

