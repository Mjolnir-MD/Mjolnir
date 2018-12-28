#ifndef MJOLNIR_READ_EXTERNAL_FORCEFIELD_HPP
#define MJOLNIR_READ_EXTERNAL_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_external_interaction.hpp>
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
ExternalForceField<traitsT>
read_external_forcefield(toml::array interactions,
                         const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_external_forcefield(), 0);
    MJOLNIR_LOG_INFO(interactions.size(),
                     " kinds of external interactions are found.");

    ExternalForceField<traitsT> eff;
    for(const auto& interaction : toml::get<std::vector<toml::table>>(interactions))
    {
        if(interaction.count("file_name") == 1)
        {
            MJOLNIR_SCOPE(interaction.count("file_name") == 1, 1);

            const auto file_name = get_toml_value<std::string>(
                        interaction, "file_name", "[[external]]");
            MJOLNIR_LOG_INFO("file_name = ", file_name);

            if(interaction.size() != 1)
            {
                MJOLNIR_LOG_WARN(
                    "[[forcefields.external]] has `file_name` and other keys.");
                MJOLNIR_LOG_WARN(
                    "When `file_name` is provided, other values are ignored "
                    "because those are read from the specified file (",
                    file_name, ").");
            }

            MJOLNIR_LOG_NOTICE("external forcefield is defined in `",
                               input_path, file_name, "`.");
            const auto forcefield_file = toml::parse(input_path + file_name);
            if(forcefield_file.count("forcefields") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[external] should be provided as a root object of file ",
                    file_name, ". but [[forcefields]] table found");
                std::exit(1);
            }
            if(forcefield_file.count("external") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[external] should be provided as a root object of file ",
                    file_name, ". but [external] key found");

                if(forcefield_file.at("external").type() != toml::value_t::Table)
                {
                    MJOLNIR_LOG_ERROR("type of `external` is different from "
                                      "toml::Table in file (", file_name, ").");
                    MJOLNIR_LOG_ERROR("note: [[...]] means Array-of-Tables. "
                                      "please take care.");
                    std::exit(1);
                }

                eff.emplace(read_external_interaction<traitsT>(
                    get_toml_value<toml::Table>(
                        forcefield_file, "external", file_name)));
            }
            else
            {
                eff.emplace(read_external_interaction<traitsT>(forcefield_file));
            }
        }
        else
        {
            eff.emplace(read_external_interaction<traitsT>(interaction));
        }
    }
    return eff;
}

} // mjolnir
#endif// MJOLNIR_READ_EXTERNAL_FORCEFIELD_HPP
