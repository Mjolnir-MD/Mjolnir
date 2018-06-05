#ifndef JARNGREIPR_READ_CHAIN_IDS
#define JARNGREIPR_READ_CHAIN_IDS
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>

namespace jarngreipr
{

template<typename charT, typename traits, typename alloc>
std::vector<std::basic_string<charT, traits, alloc>>
split_string(const std::basic_string<charT, traits, alloc>& str,
             const charT delim)
{
    std::vector<std::basic_string<charT, traits, alloc>> strs;
    std::basic_string<charT, traits, alloc> elem;
    std::basic_stringstream<charT> bss(str);
    while(std::getline(bss, elem, delim))
    {
        strs.push_back(elem);
    }
    return strs;
}

template<typename charT, typename traits, typename alloc>
std::vector<charT>
read_chain_ids(const std::basic_string<charT, traits, alloc>& key)
{
    if(key.size() == 1)
    {
        if(not std::isupper(key.front()))
        {
            throw std::runtime_error("jarngreipr::split_chain_ids: "
                    "chain ID should be specified in upper case: " + key);
        }
        return std::vector<char>(1, key.front());
    }

    std::vector<char> ids;
    for(auto elem : split_string(key, '&'))
    {
        if(elem.size() == 3 && elem.at(1) == '-')
        {
            if(not std::isupper(elem.front()) || not std::isupper(elem.back()))
            {
                throw std::runtime_error("jarngreipr::split_chain_ids: "
                        "chain ID should be specified in upper case: " + elem);
            }

            for(char c = elem.front(), e = elem.back(); c <= e; ++c)
            {
                ids.push_back(c);
            }
        }
        else if(elem.size() == 1)
        {
            if(not std::isupper(elem.front()))
            {
                throw std::runtime_error("jarngreipr::split_chain_ids: "
                        "chain ID should be specified in upper case: " + key);
            }
            ids.push_back(elem.front());
        }
        else
        {
            throw std::runtime_error("jarngreipr::split_chain_ids: "
                "invalid chain ID input (should be `X` or `X-Y` or `X&Y`): " +
                elem);
        }
    }
    return ids;
}

} // jarngreipr
#endif //JARNGREIPR_READ_CHAIN_IDS
