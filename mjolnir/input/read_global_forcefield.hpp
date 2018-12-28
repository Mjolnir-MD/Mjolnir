#ifndef MJOLNIR_READ_GLOBAL_FORCEFIELD_HPP
#define MJOLNIR_READ_GLOBAL_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_global_interaction.hpp>
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
GlobalForceField<traitsT>
read_global_forcefield(std::vector<toml::Table> interactions,
                       const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_forcefield(), 0);
    MJOLNIR_LOG_INFO(interactions.size(),
                     " kinds of local interactions are found.");

    GlobalForceField<traitsT> gff;
    for(const auto& interaction : interactions)
    {
        if(interaction.count("file_name") == 1)
        {
            MJOLNIR_SCOPE(interaction.count("file_name") == 1, 1);

            const std::string file_name = get_toml_value<std::string>(
                    interaction, "file_name", "[[global]]");
            MJOLNIR_LOG_INFO("file_name = ", file_name);

            if(interaction.size() != 1)
            {
                MJOLNIR_LOG_WARN(
                    "[[forcefields.global]] has `file_name` and other keys.");
                MJOLNIR_LOG_WARN(
                    "When `file_name` is provided, other values are ignored "
                    "because those are read from the specified file (",
                    file_name, ").");
            }

            MJOLNIR_LOG_NOTICE("global forcefield is defined in `",
                               input_path, file_name, "`.");
            const auto forcefield_file = toml::parse(input_path + file_name);
            if(forcefield_file.count("forcefields") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[global] should be provided as a root object of file ",
                    file_name, ". but [[forcefields]] table found");
                std::exit(1);
            }
            if(forcefield_file.count("global") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[global] should be provided as a root object of file ",
                    file_name, ". but [global] key found");

                if(forcefield_file.at("global").type() != toml::value_t::Table)
                {
                    MJOLNIR_LOG_ERROR("type of `global` is different from "
                                      "toml::Table in file (", file_name, ").");
                    MJOLNIR_LOG_ERROR("note: [[...]] means Array-of-Tables. "
                                      "please take care.");
                    std::exit(1);
                }
                gff.emplace(read_global_interaction<traitsT>(
                    get_toml_value<toml::Table>(
                        forcefield_file, "global", file_name)));
            }
            else
            {
                gff.emplace(read_global_interaction<traitsT>(forcefield_file));
            }
        }
        else
        {
            gff.emplace(read_global_interaction<traitsT>(interaction));
        }
    }
    return gff;
}


} // mjolnir
#endif// MJOLNIR_READ_GLOBAL_FORCEFIELD_HPP
