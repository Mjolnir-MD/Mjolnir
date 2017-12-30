#ifndef JARNGREIPR_UTIL_SPLIT_CHAIN_IDS
#define JARNGREIPR_UTIL_SPLIT_CHAIN_IDS
#include <jarngreipr/util/split_string.hpp>
#include <stdexcept>
#include <vector>
#include <cctype>

namespace mjolnir
{

template<typename charT, typename traits, typename alloc>
std::vector<charT>
split_chain_ids(const std::basic_string<charT, traits, alloc>& key)
{
    if(key.size() == 1)
    {
        if(not std::isupper(key.front()))
        {
            throw std::runtime_error("jarngreipr::split_chain_ids: "
                    "chain ID should be specified in upper case.");
        }
        return std::vector<char>{key.front()};
    }

    std::vector<char> ids;
    for(auto&& elem : split_string(key, '&'))
    {
        if(elem.size() == 3 && elem.at(1) == '-')
        {
            for(char c = std::toupper(elem.front()),
                    e = std::toupper(elem.back()); c <= e; ++c)
            {
                ids.push_back(c);
            }
        }
        else if(elem.size() == 1)
        {
            ids.push_back(std::toupper(elem.front()));
        }
        else
        {
            throw std::runtime_error("jarngreipr::split_chain_ids: "
                    "invalid chain ID input: " + elem);
        }
    }
    return ids;
}

} // mjolnir
#endif //JARNGREIPR_UTIL_SPLIT_CHAIN_IDS
