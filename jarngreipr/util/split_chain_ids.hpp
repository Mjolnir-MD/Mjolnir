#ifndef JARNGREIPR_UTIL_SPLIT_CHAIN_IDS
#define JARNGREIPR_UTIL_SPLIT_CHAIN_IDS
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cctype>

namespace mjolnir
{

inline std::vector<char> split_chain_ids(const std::string& key)
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

    if(not (key.size() == 3 && (key.at(1) == '-' || key.at(1) == '&') &&
            std::isupper(key.front()) && std::isupper(key.back())))
    {
        throw std::runtime_error("jarngreipr::split_chain_ids: "
                "chain ID must be upper case and supplied in this way: "
                "'A', 'A-C', or 'A&D'");
    }

    if(key.at(1) == '&')
    {
        // "A&D" -> {A, D}
        return std::vector<char>{key.front(), key.back()};
    }

    // "A-D" -> {A, B, C, D}
    std::vector<char> ids;
    for(char c = key.front(); c <= key.back(); ++c)
    {
        ids.push_back(c);
    }
    return ids;
}

} // mjolnir
#endif //JARNGREIPR_UTIL_SPLIT_CHAIN_IDS
