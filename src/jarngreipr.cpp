#include <extlib/toml/toml.hpp>
#include <jarngreipr/io/PDBReader.hpp>
#include <jarngreipr/io/PDBWriter.hpp>
#include <jarngreipr/io/write_as_xyz.hpp>
#include <jarngreipr/io/write_as_xyz.hpp>
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/model/make_coarse_grained.hpp>
#include <jarngreipr/forcefield/ForceFieldGenerator.hpp>
#include <jarngreipr/forcefield/ClementiGo.hpp>
// #include <jarngreipr/io/read_parameter_table.hpp>
// #include <jarngreipr/model/ClementiGo.hpp>
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <algorithm>
#include <map>

using mjolnir::operator"" _str;

namespace mjolnir
{

std::vector<char> split_chain_ids(const std::string& key)
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

template<typename coordT>
std::vector<std::unordered_map<char, PDBChain<coordT>>>
read_reference_structures(const std::vector<toml::Table>& structures)
{
    std::vector<std::unordered_map<char, PDBChain<coordT>>> tables;
    for(const auto& conf : structures)
    {
        std::unordered_map<char, PDBChain<coordT>> table;
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
                mjolnir::PDBReader<coordT> pdb_reader(filename);
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

template<typename coordT>
std::vector<std::unordered_map<char, PDBChain<coordT>>>
read_initial_structures(const std::vector<toml::Table>& structures)
{
    std::vector<std::unordered_map<char, PDBChain<coordT>>> tables;
    for(const auto& conf : structures)
    {
        std::unordered_map<char, PDBChain<coordT>> table;
        for(const auto& kvp : conf)
        {
            std::vector<char> chIDs = split_chain_ids(kvp.first);
            const toml::Table val = toml::get<toml::Table>(kvp.second);

            std::string filename;
            try
            {
                filename = toml::get<std::string>(val.at("initial"));
            }
            catch(const std::exception& except)
            {
                // do nothing. initial configuration is optional.
            }
            try
            {
                filename = toml::get<std::string>(val.at("file"));
            }
            catch(const std::exception& except)
            {
                throw std::runtime_error("jarngreipr::read_initial_structures: "
                        "neither file nor initial exists for chain " + kvp.first);
            }

            if(filename.substr(filename.size() - 4) == ".pdb")
            {
                mjolnir::PDBReader<coordT> pdb_reader(filename);
                if(not pdb_reader.is_good())
                {
                    throw std::runtime_error(
                            "jarngreipr::read_reference_structures: "
                            "file open error: " + filename);
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

inline std::vector<std::unordered_map<char, std::string>>
read_coarse_grained_models(const std::vector<toml::Table>& structures)
{
    std::vector<std::unordered_map<char, std::string>> tables;
    for(const auto& conf : structures)
    {
        std::unordered_map<char, std::string> table;
        for(const auto& kvp : conf)
        {
            const std::vector<char> chIDs = split_chain_ids(kvp.first);
            const toml::Table val = toml::get<toml::Table>(kvp.second);

            std::string model_name;
            try
            {
                model_name = toml::get<std::string>(val.at("model"));
            }
            catch(const std::exception& except)
            {
                throw std::runtime_error("jarngreipr::read_coarse_grained model: "
                        "model is not specified for chain " + kvp.first);
            }

            for(char id : chIDs)
            {
                table[id] = model_name;
            }
        }
        tables.push_back(std::move(table));
    }
    return tables;
}

template<typename coordT>
std::vector<std::unordered_map<char, CGChain<coordT>>>
apply_coarse_grained_models(
        const std::vector<std::unordered_map<char, PDBChain<coordT>>>& pdbs,
        const std::vector<std::unordered_map<char, std::string>>& models)
{
    assert(pdbs.size() == models.size());
    std::vector<std::unordered_map<char, CGChain<coordT>>> cgss;
    for(std::size_t i=0; i<pdbs.size(); ++i)
    {
        std::unordered_map<char, CGChain<coordT>> cgs;
        const auto& pdb   = pdbs.at(i);
        const auto& model = models.at(i);

        for(const auto& kv : pdb)
        {
            const auto chid = kv.first;
            const auto&  ch = kv.second;
            cgs[chid] = make_coarse_grained(ch, model.at(chid));
        }
        cgss.push_back(std::move(cgs));
    }
    return cgss;
}

} // mjolnir

int main(int argc, char **argv)
{
    typedef mjolnir::Vector<double, 3>    coord_type;
    typedef mjolnir::PDBChain<coord_type> pdb_chain_type;
    typedef mjolnir::CGChain<coord_type>  cg_chain_type;

    if(argc != 2)
    {
        std::cerr << "Usage: $ jarngreipr [file.toml]\n";
        std::cerr << "    see input/example.toml for the format" << std::endl;
        return 1;
    }

    const auto input_data = toml::parse(std::string(argv[1]));

    /* prepairing parameters */
    const auto& general  = input_data.at("general").cast<toml::value_t::Table>();
    const auto file_name = toml::get<std::string>(general.at("file_name"));
    const auto seed      = toml::get<std::int64_t>(general.at("seed"));
//     const auto parameter_table = mjolnir::read_parameters(
//             toml::get<std::string>(general.at("parameter_path")));

    /* generating coarse-grained structures */
    const auto structure_config =
            toml::get<std::vector<toml::Table>>(input_data.at("structures"));

    // vector<map<char, PDBChain<coord>>>
    // in most cases, the size of vector is one.
    const auto ref_pdbs = mjolnir::read_reference_structures<coord_type>(structure_config);
    const auto ini_pdbs = mjolnir::read_initial_structures  <coord_type>(structure_config);
    const auto models   = mjolnir::read_coarse_grained_models(structure_config);
    // Coarse Grained Structure -> cgs
    const auto ref_cgss = mjolnir::apply_coarse_grained_models(ref_pdbs, models);
    const auto ini_cgss = mjolnir::apply_coarse_grained_models(ini_pdbs, models);

    const std::unique_ptr<mjolnir::IntraChainForceFieldGenerator<coord_type>>
        ffgen = mjolnir::make_unique<mjolnir::ClementiGo<coord_type>>();

    std::size_t offset(0);
    for(const auto& chains : ref_cgss)
    {
        for(const auto& idch : chains)
        {
            const auto& ch = idch.second;
            const auto connect = ffgen->generate(std::cout, ch, offset);
            offset += ch.size();
        }
    }

    return 0;
}
