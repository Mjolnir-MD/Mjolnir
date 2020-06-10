#ifndef MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#define MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/format_nth.hpp>
#include <mjolnir/input/read_local_interaction.hpp>
#include <mjolnir/input/read_global_interaction.hpp>
#include <mjolnir/input/read_external_interaction.hpp>
#include <mjolnir/input/read_table_from_file.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/input/utility.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<ForceFieldBase<traitsT>>
read_forcefield(const toml::value& root, const std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto ff = read_table_from_file(toml::find(root, "forcefields").at(N),
                                         "forcefields");
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
        const auto& locals = toml::find(ff, "local").as_array();
        MJOLNIR_LOG_NOTICE("LocalForceField (x", locals.size(), ") found");

        for(std::size_t i=0; i<locals.size(); ++i)
        {
            MJOLNIR_LOG_NOTICE("reading ", format_nth(i), " [[forcefields.local]]");
            loc.emplace(read_local_interaction<traitsT>(
                read_table_from_file(locals.at(i), "local")));
        }
    }
    if(ff.as_table().count("global") != 0)
    {
        const auto& globals = toml::find(ff, "global").as_array();
        MJOLNIR_LOG_NOTICE("GlobalForceField (x", globals.size(), ") found");
        for(std::size_t i=0; i<globals.size(); ++i)
        {
            MJOLNIR_LOG_NOTICE("reading ", format_nth(i), " [[forcefields.global]]");
            glo.emplace(read_global_interaction<traitsT>(
                read_table_from_file(globals.at(i), "global")));
        }
    }
    if(ff.as_table().count("external") != 0)
    {
        const auto& externals = toml::find(ff, "external").as_array();
        MJOLNIR_LOG_NOTICE("ExternalForceField (x", externals.size(), ") found");
        for(std::size_t i=0; i<externals.size(); ++i)
        {
            MJOLNIR_LOG_NOTICE("reading ", format_nth(i), " [[forcefields.external]]");
            ext.emplace(read_external_interaction<traitsT>(
                read_table_from_file(externals.at(i), "external")));
        }
    }
    return make_unique<ForceField<traitsT>>(std::move(loc), std::move(glo), std::move(ext));
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_forcefield(const toml::value&, const std::size_t);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_forcefield(const toml::value&, const std::size_t);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_forcefield(const toml::value&, const std::size_t);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_forcefield(const toml::value&, const std::size_t);
#endif

} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD_HPP
