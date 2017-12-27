#ifndef JARNGREPIR_READ_COARSE_GRAINED_MODELS_HPP
#define JARNGREPIR_READ_COARSE_GRAINED_MODELS_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/get_toml_value.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <jarngreipr/io/PDBReader.hpp>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>

namespace mjolnir
{

inline std::vector<std::map<char, std::string>>
read_coarse_grained_models(const std::vector<toml::Table>& structures)
{
    std::vector<std::map<char, std::string>> tables;
    for(const auto& conf : structures)
    {
        std::map<char, std::string> table;
        for(const auto& kvp : conf)
        {
            const std::vector<char> chIDs = split_chain_ids(kvp.first);
            const toml::Table val = toml::get<toml::Table>(kvp.second);

            std::string model_name;
            model_name = toml::get<std::string>(
                    toml_value_at(val, "model", "[[structures]]"));
            for(char id : chIDs)
            {
                table[id] = model_name;
            }
        }
        tables.push_back(std::move(table));
    }
    return tables;
}

} // mjolnir
#endif //JARNGREPIR_READ_REFERENCE_STRUCTURES_HPP
