#ifndef JARNGREIPR_STRING
#define JARNGREIPR_STRING
#include <iterator>
#include <string>
#include <algorithm>
#include <cctype>

namespace mjolnir
{

template<typename charT, typename traitsT, typename allocT>
std::basic_string<charT, traitsT, allocT>
remove_whitespace(const std::basic_string<charT, traitsT, allocT>& str)
{
    std::basic_string<charT, traitsT, allocT> retval; retval.reserve(str.size());
    std::copy_if(str.cbegin(), str.cend(), std::back_inserter(retval),
                 [](const charT c){return !(std::isspace(c));});
    return retval;
}

// this is re-definition of C++14 std::literals::string_literals::operator""s.
inline std::string operator"" _str(const char* str, std::size_t len)
{
    return std::string{str, len};
}

} // mjolnir
#endif /*JARNGREIPR_STRING*/
