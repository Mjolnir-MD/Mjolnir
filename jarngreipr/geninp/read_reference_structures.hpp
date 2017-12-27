#ifndef JARNGREPIR_READ_REFERENCE_STRUCTURES_HPP
#define JARNGREPIR_READ_REFERENCE_STRUCTURES_HPP
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

template<typename coordT>
std::vector<std::map<char, PDBChain<coordT>>>
read_reference_structures(const std::vector<toml::Table>& structures)
{
    std::vector<std::map<char, PDBChain<coordT>>> tables;
    for(const auto& conf : structures)
    {
        std::map<char, PDBChain<coordT>> table;
        for(const auto& kvp : conf)
        {
            std::vector<char> chIDs = split_chain_ids(kvp.first);
            const toml::Table val = toml::get<toml::Table>(kvp.second);

            std::string filename;
            try
            {
                filename = toml::get<std::string>(val.at("file"));
            }
            catch(const std::exception& except)
            {
                throw std::runtime_error("jarngreipr::read_reference_structures: "
                        "file is not specified for chain " + kvp.first);
            }

            if(filename.substr(filename.size() - 4) == ".pdb")
            {
                PDBReader<coordT> pdb_reader(filename);
                if(not pdb_reader.is_good())
                {
                    throw std::runtime_error(
                            "jarngreipr::read_reference_structures: "
                            "file open error: filename = " + filename);
                }
                while(not pdb_reader.is_eof())
                {
                    const auto chain = pdb_reader.read_next_chain();
                    const char chain_id = chain.chain_id();
                    const auto found = std::find(
                            chIDs.begin(), chIDs.end(), chain_id);
                    if(found != chIDs.end())
                    {
                        table[*found] = chain;
                        chIDs.erase(found);
                    }
                }
            }
            else
            {
                throw std::runtime_error("jarngreipr::read_reference_structures: "
                        "unrecognizable file: " + filename);
            }

            if(not chIDs.empty())
            {
                std::string mes("jarngreipr::read_reference_structures: "
                                "missing chains in : ");
                mes += filename;
                mes += ", ID = ";
                mes += std::string(chIDs.begin(), chIDs.end());
                throw std::runtime_error(mes);
            }
        }
        tables.push_back(std::move(table));
    }
    return tables;
}
} // mjolnir
#endif // JARNGREPIR_READ_REFERENCE_STRUCTURES_HPP
