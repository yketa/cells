/*
Provides assert(expr) similar to https://en.cppreference.com/w/cpp/error/assert
but raises an std::runtime_error rather than aborting.
*/

#ifndef ASSERT_HPP
#define ASSERT_HPP

#include <stdexcept>
#include <string>

extern void __assertion_error(
    std::string const& assertion,   // asserted expression
    std::string const& func,        // function which asserted it
    std::string const& file,        // file in which this assertion occurs
    unsigned int const& line);      // line at which this assertion occurs

#endif

#undef assert
#define __FUNC __extension__ __PRETTY_FUNCTION__
#define assert(expr)                                                \
    (static_cast <bool> (expr)                                      \
        ? void (0)                                                  \
        : __assertion_error (#expr, __FUNC, __FILE__, __LINE__))

