#ifndef MJOLNIR_INPUT_READ_GLOBAL_FORCEFIELD_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_global_interaction.hpp>
#include <mjolnir/input/read_path.hpp>

namespace mjolnir
{

template<typename traitsT>
GlobalForceField<traitsT>
read_global_forcefield(const toml::array& interactions, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_INFO(interactions.size(), " global interactions are found.");

    GlobalForceField<traitsT> gff;
    for(const auto& interaction : interactions)
    {
        if(interaction.as_table().count("file_name") == 1)
        {
            MJOLNIR_LOG_SCOPE(if(interaction.as_table().count("file_name") == 1));

            const auto file_name = toml::find<std::string>(interaction, "file_name");
            if(interaction.as_table().size() != 1)
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

            const auto ff_file = toml::parse(input_path + file_name);
            if(ff_file.as_table().count("forcefield") == 1)
            {
                const auto& ff_tab = toml::find(ff_file, "forcefield");
                if(ff_tab.as_table().count("global") == 1)
                {
                    gff.emplace(read_global_interaction<traitsT>(
                            toml::find(ff_tab, "global")));
                }
            }
            throw_exception<std::runtime_error>("[error] "
                "mjolnir::read_global_forcefield: [forcefield.global] table"
                " should be provided in the file\n --> ", input_path, file_name,
                ".");
        }
        else
        {
            gff.emplace(read_global_interaction<traitsT>(interaction));
        }
    }
    return gff;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template GlobalForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_global_forcefield(const toml::array& interactions, const std::string& input_path);
extern template GlobalForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_global_forcefield(const toml::array& interactions, const std::string& input_path);
extern template GlobalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_global_forcefield(const toml::array& interactions, const std::string& input_path);
extern template GlobalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_global_forcefield(const toml::array& interactions, const std::string& input_path);
#endif

} // mjolnir
#endif// MJOLNIR_READ_GLOBAL_FORCEFIELD_HPP
