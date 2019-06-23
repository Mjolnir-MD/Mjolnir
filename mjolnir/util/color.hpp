#ifndef MJOLNIR_UTILITY_COLOR_HPP
#define MJOLNIR_UTILITY_COLOR_HPP
#include <iostream>
#include <cstdio> // ::fileno (on *nix)

#if defined(__linux__) || defined(__unix__) || defined(__APPLE__)
#include <unistd.h> // ::isatty
#endif

// This utility manipulators output some ANSI escape codes.
// On Windows (not supported by Mjolnir currently), it does not check
// the stream is connected to a tty.

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

template<typename charT, typename traits>
std::basic_ostream<charT, traits>& red(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[31m";
    }
    return os;
}
template<typename charT, typename traits>
std::basic_ostream<charT, traits>& green(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[32m";
    }
    return os;
}
template<typename charT, typename traits>
std::basic_ostream<charT, traits>& yellow(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[33m";
    }
    return os;
}
template<typename charT, typename traits>
std::basic_ostream<charT, traits>& blue(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[34m";
    }
    return os;
}
template<typename charT, typename traits>
std::basic_ostream<charT, traits>& magenta(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[35m";
    }
    return os;
}
template<typename charT, typename traits>
std::basic_ostream<charT, traits>& cyan(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[36m";
    }
    return os;
}
template<typename charT, typename traits>
std::basic_ostream<charT, traits>& white(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[37m";
    }
    return os;
}
template<typename charT, typename traits>
std::basic_ostream<charT, traits>& nocolor(std::basic_ostream<charT, traits>& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[0m";
    }
    return os;
}

template<typename charT, typename traits>
struct basic_lock_nocolor
{
    basic_lock_nocolor(std::basic_ostream<charT, traits>& os) noexcept
        : os_(std::addressof(os))
    {}
    ~basic_lock_nocolor() noexcept
    {
        if(os_)
        {
            *os_ << nocolor;
        }
    }
  private:
    std::basic_ostream<charT, traits>* os_;
};

template<typename charT, typename traits>
basic_lock_nocolor<charT, traits>
lock_nocolor(std::basic_ostream<charT, traits>& os) noexcept
{
    return basic_lock_nocolor<charT, traits>(os);
}

} // io
} // mjolnir
#endif// MJOLNIR_UTILITY_COLOR_HPP
