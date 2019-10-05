#ifndef MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#define MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_interaction.hpp>
#include <mjolnir/input/read_global_interaction.hpp>
#include <mjolnir/input/read_external_interaction.hpp>
#include <mjolnir/input/read_table_from_file.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/input/utility.hpp>

namespace mjolnir
{

template<typename traitsT>
ForceField<traitsT>
read_forcefield(const toml::value& root, const std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto input_path = read_input_path(root);
    const auto ff = read_table_from_file(root, "forcefields", N, input_path);
    check_keys_available(ff, {"local"_s, "global"_s, "external"_s, "name"_s});

    if(ff.as_table().count("name") == 1)
    {
        MJOLNIR_LOG_NOTICE("forcefield \"", toml::find<std::string>(ff, "name"),
                           "\" found");
    }

    LocalForceField<traitsT>    loc;
    GlobalForceField<traitsT>   glo;
    ExternalForceField<traitsT> ext;

    if(ff.as_table().count("local") != 0)
    {
        const auto num_loc = toml::find(ff, "local").as_array().size();
        MJOLNIR_LOG_INFO("LocalForceField (x", num_loc, ") found");
        for(std::size_t i=0; i<num_loc; ++i)
        {
            loc.emplace(read_local_interaction<traitsT>(
                read_table_from_file(ff, "local", i, input_path)));
        }
    }
    if(ff.as_table().count("global") != 0)
    {
        const auto num_glo = toml::find(ff, "global").as_array().size();
        MJOLNIR_LOG_INFO("GlobalForceField (x", num_glo, ") found");
        for(std::size_t i=0; i<num_glo; ++i)
        {
            glo.emplace(read_global_interaction<traitsT>(
                read_table_from_file(ff, "global", i, input_path)));
        }
    }
    if(ff.as_table().count("external") != 0)
    {
        const auto num_ext = toml::find(ff, "external").as_array().size();
        MJOLNIR_LOG_INFO("ExternalForceField (x", num_ext, ") found");
        for(std::size_t i=0; i<num_ext; ++i)
        {
            ext.emplace(read_external_interaction<traitsT>(
                read_table_from_file(ff, "external", i, input_path)));
        }
    }
    return ForceField<traitsT>(std::move(loc), std::move(glo), std::move(ext));
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template ForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_forcefield(const toml::value& root, const std::size_t N);
extern template ForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_forcefield(const toml::value& root, const std::size_t N);
extern template ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, const std::size_t N);
extern template ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, const std::size_t N);
#endif

} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD_HPP
