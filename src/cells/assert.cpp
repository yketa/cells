#include <stdexcept>
#include <string>

#include "assert.hpp"

extern void __assertion_error(
    std::string const& assertion,
    std::string const& func,
    std::string const& file,
    unsigned int const& line) {
    // error message
    std::string msg =
        file + ":" + std::to_string(line) + ": " + func
            + ": Assertion `" + assertion + "' failed.";
    // throw error
    throw std::runtime_error(msg);
}

