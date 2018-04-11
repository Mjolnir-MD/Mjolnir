#ifndef MJOLNIR_STRING_HPP
#define MJOLNIR_STRING_HPP
#include <string>

namespace mjolnir
{

inline std::string operator"" _str(const char* str, std::size_t len)
{
    return std::string{str, len};
}

} // mjolnir
#endif// MJOLNIR_STRING_HPP
