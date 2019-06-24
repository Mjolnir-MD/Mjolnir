#ifndef MJOLNIR_UTILITY_IO_HPP
#define MJOLNIR_UTILITY_IO_HPP
#include <memory>
#include <iostream>
#include <cstdio> // ::fileno (on *nix)
#if defined(__linux__) || defined(__unix__) || defined(__APPLE__)
#include <unistd.h> // ::isatty
#endif

// A set of I/O utility functions

namespace mjolnir
{
namespace io
{
namespace detail
{

template<typename charT, typename traits>
bool isatty(const std::basic_ostream<charT, traits>& os) noexcept
{
#if defined(__linux__) || defined(__unix__) || defined(__APPLE__)
    if(std::addressof(os) == std::addressof(std::cout))
    {
        return ::isatty(::fileno(stdout)) != 0;
    }
    else if(std::addressof(os) == std::addressof(std::cerr))
    {
        return ::isatty(::fileno(stderr)) != 0;
    }
    return false;
#elif defined(_WIN32)
    // TODO? When I bought windows machine...
    return true;
#else
#  error "unknown platform"
#endif
}

} // detail
} // io
} // mjolnir
#endif //MJOLNIR_UTILITY_IO_HPP
