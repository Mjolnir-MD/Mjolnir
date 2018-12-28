#ifndef MJOLNIR_READ_FORCEFIELD_HPP
#define MJOLNIR_READ_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_forcefield.hpp>
#include <mjolnir/input/read_global_forcefield.hpp>
#include <mjolnir/input/read_external_forcefield.hpp>
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
ForceField<traitsT>
read_forcefield_from_table(const toml::value& ff, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_forcefield_from_table(), 0);

    toml::array fflocal, ffglobal, ffexternal;
    if(toml::get<toml::table>(ff).count("local") == 1)
    {
        MJOLNIR_LOG_INFO("LocalForceField found");
        fflocal = toml::find<toml::array>(ff, "local");
    }
    if(toml::get<toml::table>(ff).count("global") == 1)
    {
        MJOLNIR_LOG_INFO("GlobalForceField found");
        ffglobal = toml::find<toml::array>(ff, "global");
    }
    if(toml::get<toml::table>(ff).count("external") == 1)
    {
        MJOLNIR_LOG_INFO("ExternalForceField found");
        ffexternal = toml::find<toml::array>(ff, "external");
    }

    for(const auto kv: toml::get<toml::table>(ff))
    {
        if(kv.first != "local" && kv.first != "global" && kv.first != "external")
        {
            MJOLNIR_LOG_WARN("unknown key `", kv.first, "` appeared in "
                             "[[forcefields]] table. It will be ignored.");
        }
    }
    return ForceField<traitsT>(
        read_local_forcefield<traitsT>(std::move(fflocal), input_path),
        read_global_forcefield<traitsT>(std::move(ffglobal), input_path),
        read_external_forcefield<traitsT>(std::move(ffexternal), input_path));
}


template<typename traitsT>
ForceField<traitsT>
read_forcefield(const toml::table& root, std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_forcefield(), 0);
    if(N != 0) {MJOLNIR_LOG_NOTICE("reading ", N, "-th [[forcefield]].");}
    else       {MJOLNIR_LOG_NOTICE("reading [[forcefield]].");}

    const auto& ffs = toml::find<toml::array>(root, "forcefields", "<root>");
    if(ffs.size() <= N)
    {
        throw_exception<std::out_of_range>("[error] mjolnir::read_forcefield: "
            "no enough forcefield definitions: ", N, " is required, but only ",
            ffs.size(), " defined.");
    }
    MJOLNIR_LOG_INFO(ffs.size(), " forcefields are provided");
    MJOLNIR_LOG_INFO("using ", N, "-th forcefield");

    const auto input_path = read_input_path(root);
    const auto& ff = toml::get<toml::table>(ffs.at(N));
    if(ff.count("file_name") == 1)
    {
        MJOLNIR_SCOPE(ff.count("file_name") == 1, 1);

        const auto file_name = toml::find<std::string>(ff, "file_name");
        MJOLNIR_LOG_NOTICE("forcefield is defined in ", input_path, file_name);
        if(ff.size() != 1)
        {
            MJOLNIR_LOG_WARN("[[forcefields]] has \"file_name\" and other keys.");
            MJOLNIR_LOG_WARN("When \"file_name\" is provided, other values are "
                             "ignored because those are read from the specified"
                             " file (", input_path, file_name, ").");
        }
        const auto ff_file = toml::parse(input_path + file_name);
        if(ff_file.count("forcefields") != 1)
        {
            throw_exception<std::out_of_range>("[error] mjolnir::read_system: "
                "table [system] not found in toml file\n --> ",
                input_path, file_name, "\n | the file should define [forcefield] "
                "table and define values in it.");
        }
        return read_forcefield_from_table<traitsT>(
                toml::find(ff_file, "forcefields"), input_path);
    }
    // all-in-one file
    return read_forcefield_from_table<traitsT>(ffs.at(N), input_path);
}

} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD_HPP
