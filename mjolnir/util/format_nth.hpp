#ifndef MJOLNIR_UTIL_FORMAT_NTH_HPP
#define MJOLNIR_UTIL_FORMAT_NTH_HPP
#include <mjolnir/util/string.hpp>
#include <string>

namespace mjolnir
{
template<typename intT>
std::string format_nth(const intT idx)
{
    if(10 < idx && idx < 20) {return std::to_string(idx) + "th"_s;}

    const auto idx_s = std::to_string(idx);
    switch(idx_s.back())
    {
        case '1': {return idx_s + "st";}
        case '2': {return idx_s + "nd";}
        case '3': {return idx_s + "rd";}
        default:  {return idx_s + "th";}
    }
}
} // mjolnir
#endif // MJOLNIR_UTIL_FORMAT_NTH_HPP
