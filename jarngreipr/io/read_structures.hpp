#ifndef JARNGREIPR_READ_STRUCTURES_H
#define JARNGREIPR_READ_STRUCTURES_H
#include <jarngriepr/model/CGChain.hpp>
#include <jarngriepr/model/CarbonAlpha.hpp>
#include <jarngriepr/io/read_chain_ids.hpp>
#include <mjolnir/util/get_toml_values.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace jarngreipr
{

template<typename realT>
std::vector<PDBChain<realT>>
read_reference_structures(const toml::Table& data, const std::size_t N)
{
    const toml::Array& structures = mjolnir::toml_value_at(
            data, "structures", "<root>").cast<toml::value_t::Array>();
    if(structures.size() <= N)
    {
        mjolnir::throw_exception<std::out_of_range>("jarngreipr::"
            "read_reference_structures: there are no structure ", N, " in the "
            "input file ([[structures]]).");
    }

    const toml::Table& structure = structures.at(N).cast<toml::value_t::Table>();

    std::vector<PDBChain<realT>> pdbs;
    for(const auto& kv : structure)
    {
        const std::string& k = kv.first;
        const toml::Table& v = kv.second.cast<toml::value_t::Table>();
        if(v.count("file") == 0)
        {
            mjolnir::throw_exception<std::out_of_range>("jarngreipr::"
                "read_reference_structures: `file` is not specified in the "
                "input file [[structures.", k, "]].");
        }
        const std::string pdbfilename = toml::get<std::string>(v.at("file"));
        PDBReader<double> reader(pdbfilename);

        const std::vector<char> chain_ids = read_chain_ids(k);
        for(const char id : chain_ids)
        {
            pdbs.push_back(reader.read_chain(id));
        }
    }
    return pdbs;
}

template<typename realT>
std::vector<PDBChain<realT>>
read_initial_structures(const toml::Table& data)
{
    const toml::Array& structures = mjolnir::toml_value_at(
            data, "structures", "<root>").cast<toml::value_t::Array>();
    if(structures.size() <= N)
    {
        mjolnir::throw_exception<std::out_of_range>("jarngreipr::"
            "read_initial_structures: there are no structure ", N, " in the "
            "input file ([[structures]])");
    }

    const toml::Table& structure = structures.at(N).cast<toml::value_t::Table>();

    std::vector<PDBChain<realT>> pdbs;
    for(const auto& kv : structure)
    {
        const std::string& k = kv.first;
        const toml::Table& v = kv.second.cast<toml::value_t::Table>();

        std::string pdbfilename;
        if(v.count("initial") == 1)
        {
            pdbfilename = toml::get<std::string>(v.at("initial"));
        }
        else if(v.count("file") == 1)
        {
            pdbfilename = toml::get<std::string>(v.at("file"));
        }
        else
        {
            mjolnir::throw_exception<std::out_of_range>("jarngreipr::"
                "read_initial_structures: neither `file` or `initial` is "
                "specified in the input file [[structures.", k, "]]");
        }
        PDBReader<double> reader(pdbfilename);

        const std::vector<char> chain_ids = read_chain_ids(k);
        for(const char id : chain_ids)
        {
            pdbs.push_back(reader.read_chain(id));
        }
    }
    return pdbs;
}
} // jarngreipr
#endif// JARNGREIPR_READ_STRUCTURES_H
