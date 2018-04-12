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

} // mjolnir
#endif /*JARNGREIPR_STRING*/
