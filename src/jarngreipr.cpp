#include <extlib/toml/toml.hpp>
#include <jarngreipr/io/PDBReader.hpp>
#include <jarngreipr/io/PDBWriter.hpp>
#include <jarngreipr/io/write_as_xyz.hpp>
#include <jarngreipr/io/write_as_xyz.hpp>
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/model/make_coarse_grained.hpp>
// #include <jarngreipr/io/read_parameter_table.hpp>
// #include <jarngreipr/model/ClementiGo.hpp>
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <algorithm>
#include <map>

using mjolnir::operator"" _str;

std::vector<char> split_ids(const std::string& key)
{
    std::vector<char> ids;
    if(key.find('-') == std::string::npos)
    {
        ids.push_back(key.front());
    }
    else
    {
        const char first = key.front();
        const char last  = key.back();
        for(char c = first; c <= last; ++c)
        {
            ids.push_back(c);
        }
    }
    return ids;
}

int main(int argc, char **argv)
{
    typedef mjolnir::Vector<double, 3>    coord_type;
    typedef mjolnir::PDBChain<coord_type> pdb_chain_type;
    typedef mjolnir::CGChain<coord_type>  cg_chain_type;

    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " [file.toml]\n";
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
    const auto& structure_tables = toml::get<std::vector<toml::Table>>(
            input_data.at("structures"));

    // structure = listof{chain}; chain = pairof{reference, initial}
    std::vector<std::vector<std::tuple<char, cg_chain_type, cg_chain_type>>>
        structure_sets;

    std::size_t structure_index=0;
    for(const auto& structure_table : structure_tables)
    {
        std::vector<std::tuple<char, cg_chain_type, cg_chain_type>>
            structure_set;

        if(structure_table.count("index") == 0)
        {
            ++structure_index; // XXX more safer way?
        }
        else
        {
            structure_index =
                toml::get<std::size_t>(structure_table.at("index"));
        }

        for(const auto& kv : structure_table)
        {
            // value = {file = "*.pdb", initial = "*.pdb", model = "CarbonAlpha"}
            const auto& value = kv.second.cast<toml::value_t::Table>();
            for(char cid : split_ids(kv.first))
            {
                structure_set.emplace_back(
                        cid, cg_chain_type{}, cg_chain_type{});
            }
            const auto model = toml::get<std::string>(value.at("model"));

            // read pdb file
            const std::string pdb_filename =
                toml::get<std::string>(value.at("file"));

            std::string ini_filename = pdb_filename;
            try
            {
                ini_filename = toml::get<std::string>(value.at("initial"));
            }
            catch(...)
            {
                // use pdb_filename as initial configuration file.
            }
            /* read reference file */{
                mjolnir::PDBReader<coord_type> pdb_reader(pdb_filename);
                if(not pdb_reader.is_good())
                {
                    std::cerr << "jarngreipr: at main(): file open error \n"
                        << "    while reading pdb file in [[structures]] table.: "
                        << pdb_filename << std::endl;
                    return 1;
                }

                // find all chains
                while(not pdb_reader.is_eof())
                {
                    const auto chain = pdb_reader.read_next_chain();
                    const char chain_id = chain.chain_id();
                    const auto found = std::find_if(
                        structure_set.begin(), structure_set.end(),
                        [=](const std::tuple<char, cg_chain_type, cg_chain_type>& v){
                            return std::get<0>(v) == chain_id;
                        });
                    if(found != structure_set.end())
                    {
                        std::get<1>(*found) =
                            mjolnir::make_coarse_grained(chain, model);
                    }
                }
                for(const auto& item : structure_set)
                {
                    if(std::get<1>(item).empty())
                    {
                        throw std::runtime_error(
                            "jarngreipr: at main(): missing chain "_str +
                            std::string(1, std::get<0>(item)) + " in file "_str +
                            ini_filename);
                    }
                }
            }
            /* read reference file */{
                mjolnir::PDBReader<coord_type> ini_reader(ini_filename);
                if(not ini_reader.is_good())
                {
                    std::cerr << "jarngreipr: at main(): file open error \n"
                        << "    while reading pdb file in [[structures]] table.: "
                        << ini_filename << std::endl;
                    return 1;
                }

                // find all chains
                while(not ini_reader.is_eof())
                {
                    const auto chain = ini_reader.read_next_chain();
                    const char chain_id = chain.chain_id();
                    const auto found = std::find_if(
                        structure_set.begin(), structure_set.end(),
                        [=](const std::tuple<char, cg_chain_type, cg_chain_type>& v){
                            return std::get<0>(v) == chain_id;
                        });
                    if(found != structure_set.end())
                    {
                        std::get<2>(*found) =
                            mjolnir::make_coarse_grained(chain, model);
                    }
                }
                for(const auto& item : structure_set)
                {
                    if(std::get<2>(item).empty())
                    {
                        throw std::runtime_error(
                            "jarngreipr: at main(): missing chain "_str +
                            std::string(1, std::get<0>(item)) + " in file "_str +
                            ini_filename);
                    }
                }
            }
        }
        structure_sets.push_back(std::move(structure_set));
    }
    /* output coarse-grained structure */{
        for(const auto& structure_set : structure_sets)
        {
            std::ofstream refstr(file_name + "_cg_"_str +
                                 std::to_string(structure_index) + ".xyz"_str);
            std::ofstream inicon(file_name + "_init_"_str +
                                 std::to_string(structure_index) + ".xyz"_str);
            std::size_t offset = 0;
            for(const auto& structure : structure_set)
            {
                write_as_pdb(refstr, std::get<1>(structure), offset);
                write_as_pdb(inicon, std::get<2>(structure), offset);
                offset += structure.size();
            }
        }
    }






    return 0;
}
