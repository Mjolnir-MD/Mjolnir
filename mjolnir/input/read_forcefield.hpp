#ifndef MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#define MJOLNIR_INPUT_READ_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasinForceField.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasin2BasinUnit.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasin3BasinUnit.hpp>
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
read_default_forcefield(const toml::value& root, const std::size_t N)
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

template<typename traitsT>
std::unique_ptr<ForceFieldBase<traitsT>>
read_multiple_basin_forcefield(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using namespace mjolnir::literals::string_literals;
    using real_type = typename traitsT::real_type;

    // ```toml
    // [simulator]
    //
    // forcefields.type = "MultipleBasin"
    //
    // [[simulator.forcefields.units]]
    // basins   = ["apo1", "open1", "close1"]
    // dVs      = [   0.0,     1.2,      6.0]
    // delta.apo1-open1   = 6.0
    // delta.open1-close1 = 5.0
    // delta.close1-apo1  = 4.0
    //
    // [[simulator.forcefields.units]]
    // basins = ["apo2", "open2", "close2"]
    // dVs    = [   0.0,     1.2,      6.0]
    // delta.apo2-open2   = 6.0
    // delta.open2-close2 = 5.0
    // delta.close2-apo2  = 4.0
    //
    // [[forcefields]]
    // name = "open1"
    // [[forcefields.local]]
    // interaction = "BondLength"
    // # ...
    //
    // [[forcefields]]
    // name = "close1"
    // [[forcefields.local]]
    // interaction = "BondLength"
    // # ...
    //
    // [[forcefields]]
    // name = "common" # "common" forcefield will be applied in any case.
    // [[forcefields.global]]
    // interaction = "Pair"
    // # ...
    // ```

    // read all the [[forcefields]] tables and make a map name -> [[forcefields]]
    std::map<std::string, std::pair<bool, toml::value const*>> ffs;
    for(const auto& ff : root.at("forcefields").as_array())
    {
        check_keys_available(ff, {"local"_s, "global"_s, "external"_s, "name"_s});
        if(!ff.contains("name"))
        {
            MJOLNIR_LOG_ERROR("MultipleBasin requires name of forcefields.");
        }
        const auto name = toml::find<std::string>(ff, "name");
        MJOLNIR_LOG_NOTICE("forcefield \"", name, "\" found");
        ffs[name] = std::make_pair(false, std::addressof(ff));
    }

    // a helper function to read (loc, glo, ext) forcefields
    const auto read_forcefield_elements = [&](const toml::value* ffptr)
    {
        LocalForceField<traitsT>    loc;
        GlobalForceField<traitsT>   glo;
        ExternalForceField<traitsT> ext;
        if(!ffptr)
        {
            return std::make_tuple(loc, glo, ext);
        }
        const auto& ff = *ffptr;

        if(ff.contains("local"))
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
        if(ff.contains("global"))
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
        if(ff.contains("external"))
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
        return std::make_tuple(std::move(loc), std::move(glo), std::move(ext));
    };

    std::vector<std::unique_ptr<MultipleBasinUnitBase<traitsT>>> units;
    for(const auto& unit : simulator.at("forcefields").at("units").as_array())
    {
        const auto& names = toml::find<std::vector<std::string>>(unit, "basins");
        if(names.size() == 2)
        {
            const auto dVs = toml::find<std::vector<real_type>>(unit, "dVs");
            if(dVs.size() != 2)
            {
                MJOLNIR_LOG_ERROR("MultipleBasin requires exactly the same "
                                  "number of dVs as that of basins.");
                throw std::runtime_error(toml::format_error("mjolnir::"
                    "read_multiple_basin_forcefield: different number of dVs "
                    "and basins.", unit.at("basins"), "2 basins are defined.",
                    unit.at("dVs"), "2 dVs expected."));
            }
            MJOLNIR_LOG_NOTICE("dV    = ", dVs);
            const auto delta = -std::abs(toml::find<real_type>(unit, "delta"));
            MJOLNIR_LOG_NOTICE("delta = ", delta);

            if(ffs.at(names.at(0)).first)
            {
                MJOLNIR_LOG_WARN("forcefield ", names.at(0), " is used more "
                                 "than twice in MultipleBasin.");
            }
            if(ffs.at(names.at(1)).first)
            {
                MJOLNIR_LOG_WARN("forcefield ", names.at(1), " is used more "
                                 "than twice in MultipleBasin.");
            }
            ffs.at(names.at(0)).first = true;
            ffs.at(names.at(1)).first = true;

            auto ff1 = read_forcefield_elements(ffs.at(names.at(0)).second);
            auto ff2 = read_forcefield_elements(ffs.at(names.at(1)).second);

            units.push_back(make_unique<MultipleBasin2BasinUnit<traitsT>>(delta,
                    names.at(0), names.at(1), dVs.at(0), dVs.at(1),
                    std::move(ff1), std::move(ff2)));
        }
        else if(names.size() == 3)
        {
            const auto dVs = toml::find<std::array<real_type, 3>>(unit, "dVs");
            if(dVs.size() != 3)
            {
                MJOLNIR_LOG_ERROR("MultipleBasin requires exactly the same "
                                  "number of dVs as that of basins.");
                throw std::runtime_error(toml::format_error("mjolnir::"
                    "read_multiple_basin_forcefield: different number of dVs "
                    "and basins.", unit.at("basins"), "3 basins are defined.",
                    unit.at("dVs"), "3 dVs expected."));
            }
            MJOLNIR_LOG_NOTICE("dV    = ", dVs);

            const auto& delta = unit.at("delta");
            const auto delta12 = toml::find<real_type>(delta, names.at(0) + "-"_s + names.at(1));
            const auto delta23 = toml::find<real_type>(delta, names.at(1) + "-"_s + names.at(2));
            const auto delta31 = toml::find<real_type>(delta, names.at(2) + "-"_s + names.at(0));
            MJOLNIR_LOG_NOTICE("delta12 = ", delta12);
            MJOLNIR_LOG_NOTICE("delta23 = ", delta23);
            MJOLNIR_LOG_NOTICE("delta31 = ", delta31);

            if(ffs.at(names.at(0)).first)
            {
                MJOLNIR_LOG_WARN("forcefield ", names.at(0), " is used more "
                                 "than twice in MultipleBasin.");
            }
            if(ffs.at(names.at(1)).first)
            {
                MJOLNIR_LOG_WARN("forcefield ", names.at(1), " is used more "
                                 "than twice in MultipleBasin.");
            }
            if(ffs.at(names.at(2)).first)
            {
                MJOLNIR_LOG_WARN("forcefield ", names.at(2), " is used more "
                                 "than twice in MultipleBasin.");
            }
            auto ff1 = read_forcefield_elements(ffs.at(names.at(0)).second);
            auto ff2 = read_forcefield_elements(ffs.at(names.at(1)).second);
            auto ff3 = read_forcefield_elements(ffs.at(names.at(2)).second);

            units.push_back(make_unique<MultipleBasin3BasinUnit<traitsT>>(
                    names.at(0), names.at(1), names.at(2),
                    delta12, delta23, delta31, dVs.at(0), dVs.at(1), dVs.at(2),
                    std::move(ff1), std::move(ff2), std::move(ff3)));
        }
        else
        {
            MJOLNIR_LOG_ERROR("MultipleBasin with 4 or more basins is not supported now.");
            throw std::runtime_error("mjolnir::read_multiple_basin_forcefield: "
                                     "too many basins.");
        }
    }

    toml::value const* common = nullptr;
    for(const auto& ff: ffs)
    {
        if(ff.first == "common")
        {
            common = ff.second.second;
            if(ff.second.first)
            {
                MJOLNIR_LOG_WARN("\"common\" forcefield is used as a basin.");
                MJOLNIR_LOG_WARN("\"common\" forcefield is considered as a "
                    "common part and applied regardress of MultipleBasin states."
                    "It may cause double-counting.");
            }
            continue;
        }
        if(ff.second.first)
        {
            continue; // used. good.
        }
        else // non-"common" forcefield is defined but not used.
        {
            MJOLNIR_LOG_WARN("forcefield \"", ff.first, "\" is not referenced "
                "from any MultipleBasin unit. Did you forget to use this?");
        }
    }
    return make_unique<MultipleBasinForceField<traitsT>>(
            read_forcefield_elements(common), std::move(units));
}


template<typename traitsT>
std::unique_ptr<ForceFieldBase<traitsT>>
read_forcefield(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if(!simulator.contains("forcefields"))
    {
        return read_default_forcefield<traitsT>(root, 0);
    }
    else if(simulator.at("forcefields").at("type").as_string() == "MultipleBasin")
    {
        return read_multiple_basin_forcefield<traitsT>(root, simulator);
    }
    else
    {
        throw std::runtime_error(toml::format_error("mjolnir::read_forcefield: "
            "unknown forcefield type", simulator.at("forcefields").at("type"),
            "here", {"expected one of the following: ",
                "- \"MultipleBasin\": Multiple Basin forcefield.",
                "- (nothing)      : In case of normal forcefield, you don't need this field."
            }));
    }
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_default_forcefield(const toml::value&, const std::size_t);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_default_forcefield(const toml::value&, const std::size_t);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_default_forcefield(const toml::value&, const std::size_t);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_default_forcefield(const toml::value&, const std::size_t);

extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_multiple_basin_forcefield(const toml::value&, const toml::value&);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_multiple_basin_forcefield(const toml::value&, const toml::value&);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_multiple_basin_forcefield(const toml::value&, const toml::value&);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_multiple_basin_forcefield(const toml::value&, const toml::value&);

extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_forcefield(const toml::value&, const toml::value&);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_forcefield(const toml::value&, const toml::value&);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_forcefield(const toml::value&, const toml::value&);
extern template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_forcefield(const toml::value&, const toml::value&);
#endif

} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD_HPP
