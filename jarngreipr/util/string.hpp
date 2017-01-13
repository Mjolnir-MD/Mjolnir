#ifndef JARNGREIPR_STRING
#define JARNGREIPR_STRING
#include <string>
#include <cctype>

namespace jarngreipr
{

template<typename charT, typename traitsT, typename allocT>
std::basic_string<charT, traitsT, allocT>
remove_whitespace(const std::basic_string<charT, traitsT, allocT>& str)
{
    std::basic_string<charT, traitsT, allocT> retval;
    for(auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        if(!std::isspace(*iter))
            retval.push_back(*iter);
    }
    return retval;
}

} // jarngreipr
#endif /*JARNGREIPR_STRING*/
