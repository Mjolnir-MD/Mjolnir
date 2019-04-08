#ifndef MJOLNIR_INPUT_READ_LOCAL_FORCEFIELD_HPP
#define MJOLNIR_INPUT_READ_LOCAL_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_interaction.hpp>
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
LocalForceField<traitsT>
read_local_forcefield(
        const toml::array& interactions, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_INFO(interactions.size(), " local interactions are found.");

    LocalForceField<traitsT> lff;
    for(const auto& interaction : interactions)
    {
        if(toml::get<toml::table>(interaction).count("file_name") == 1)
        {
            MJOLNIR_LOG_SCOPE(if(toml::get<toml::table>(interaction).count("file_name") == 1));

            const auto file_name = toml::find<std::string>(interaction, "file_name");
            if(toml::get<toml::table>(interaction).size() != 1)
            {
                MJOLNIR_LOG_WARN(
                    "[[forcefields.local]] has `file_name` and other keys.");
                MJOLNIR_LOG_WARN(
                    "When `file_name` is provided, other values are ignored "
                    "because those are read from the specified file (",
                    input_path, file_name, ").");
            }
            MJOLNIR_LOG_NOTICE("local forcefield is defined in `",
                               input_path, file_name, "`.");

            const auto ff_file = toml::parse(input_path + file_name);
            if(ff_file.count("forcefield") == 1)
            {
                const auto& ff_tab = toml::find(ff_file, "forcefield");
                if(toml::get<toml::table>(ff_tab).count("local") == 1)
                {
                    lff.emplace(read_local_interaction<traitsT>(
                            toml::find(ff_tab, "local")));
                }
            }
            throw_exception<std::runtime_error>("[error] "
                "mjolnir::read_local_forcefield: [forcefield.local] table"
                " should be provided in the file\n --> ", input_path, file_name,
                ".");
        }
        else
        {
            lff.emplace(read_local_interaction<traitsT>(interaction));
        }
    }
    return lff;
}

} // mjolnir
#endif// MJOLNIR_READ_LOCAL_FORCEFIELD_HPP
