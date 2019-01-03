#ifndef MJOLNIR_UTIL_STRING_HPP
#define MJOLNIR_UTIL_STRING_HPP
#include <string>

namespace mjolnir
{
inline namespace literals
{
inline namespace string_literals
{

inline std::string operator"" _s(const char* str, std::size_t len)
{
    return std::string{str, len};
}

} // string_literals
} // literals
} // mjolnir
#endif// MJOLNIR_STRING_HPP
