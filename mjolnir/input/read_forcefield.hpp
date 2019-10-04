#ifndef MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#define MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_forcefield.hpp>
#include <mjolnir/input/read_global_forcefield.hpp>
#include <mjolnir/input/read_external_forcefield.hpp>
#include <mjolnir/input/read_table_from_file.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/input/utility.hpp>

namespace mjolnir
{

template<typename traitsT>
ForceField<traitsT>
read_forcefield_from_table(const toml::value& ff, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    check_keys_available(ff, {"local"_s, "global"_s, "external"_s, "name"_s});

    if(ff.as_table().count("name") == 1)
    {
        MJOLNIR_LOG_NOTICE("forcefield \"", toml::find<std::string>(ff, "name"),
                           "\" found");
    }

    toml::array fflocal, ffglobal, ffexternal;
    if(ff.as_table().count("local") == 1)
    {
        MJOLNIR_LOG_INFO("LocalForceField found");
        fflocal = toml::find<toml::array>(ff, "local");
    }
    if(ff.as_table().count("global") == 1)
    {
        MJOLNIR_LOG_INFO("GlobalForceField found");
        ffglobal = toml::find<toml::array>(ff, "global");
    }
    if(ff.as_table().count("external") == 1)
    {
        MJOLNIR_LOG_INFO("ExternalForceField found");
        ffexternal = toml::find<toml::array>(ff, "external");
    }

    return ForceField<traitsT>(
        read_local_forcefield<traitsT>(std::move(fflocal), input_path),
        read_global_forcefield<traitsT>(std::move(ffglobal), input_path),
        read_external_forcefield<traitsT>(std::move(ffexternal), input_path));
}

template<typename traitsT>
ForceField<traitsT>
read_forcefield(const toml::value& root, const std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto ff = read_table_from_file(
            root, "forcefields", N, read_input_path(root));

    check_keys_available(ff, {"local"_s, "global"_s, "external"_s, "name"_s});

    if(ff.as_table().count("name") == 1)
    {
        MJOLNIR_LOG_NOTICE("forcefield \"", toml::find<std::string>(ff, "name"),
                           "\" found");
    }

    toml::array fflocal, ffglobal, ffexternal;
    if(ff.as_table().count("local") == 1)
    {
        MJOLNIR_LOG_INFO("LocalForceField found");
        fflocal = toml::find<toml::array>(ff, "local");
    }
    if(ff.as_table().count("global") == 1)
    {
        MJOLNIR_LOG_INFO("GlobalForceField found");
        ffglobal = toml::find<toml::array>(ff, "global");
    }
    if(ff.as_table().count("external") == 1)
    {
        MJOLNIR_LOG_INFO("ExternalForceField found");
        ffexternal = toml::find<toml::array>(ff, "external");
    }

    const auto input_path = read_input_path(root);

    return ForceField<traitsT>(
        read_local_forcefield<traitsT>(std::move(fflocal), input_path),
        read_global_forcefield<traitsT>(std::move(ffglobal), input_path),
        read_external_forcefield<traitsT>(std::move(ffexternal), input_path));
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template ForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_forcefield_from_table(const toml::value& ff, const std::string& input_path);
extern template ForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_forcefield_from_table(const toml::value& ff, const std::string& input_path);
extern template ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_forcefield_from_table(const toml::value& ff, const std::string& input_path);
extern template ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_forcefield_from_table(const toml::value& ff, const std::string& input_path);

extern template ForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_forcefield(const toml::value& root, std::size_t N);
extern template ForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_forcefield(const toml::value& root, std::size_t N);
extern template ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, std::size_t N);
extern template ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, std::size_t N);
#endif

} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD_HPP
