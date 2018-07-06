#ifndef MJOLNIR_READ_FORCEFIELD
#define MJOLNIR_READ_FORCEFIELD
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/input/read_interaction.hpp>

namespace mjolnir
{

template<typename traitsT>
LocalForceField<traitsT>
read_local_forcefield(std::vector<toml::Table> interactions)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_local_forcefield(), 0);
    LocalForceField<traitsT> lff;
    for(const auto& interaction : interactions)
    {
        lff.emplace(read_local_interaction<traitsT>(interaction));
    }
    return lff;
}

template<typename traitsT>
GlobalForceField<traitsT>
read_global_forcefield(std::vector<toml::Table> interactions)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_forcefield(), 0);
    GlobalForceField<traitsT> gff;
    for(const auto& interaction : interactions)
    {
        gff.emplace(read_global_interaction<traitsT>(interaction));
    }
    return gff;
}

template<typename traitsT>
ExternalForceField<traitsT>
read_external_forcefield(std::vector<toml::Table> interactions)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_external_forcefield(), 0);
    ExternalForceField<traitsT> eff;
    for(const auto& interaction : interactions)
    {
        eff.emplace(read_external_interaction<traitsT>(interaction));
    }
    return eff;
}

template<typename traitsT>
ForceField<traitsT>
read_forcefield_from_table(const toml::Table& ff)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_forcefield_from_table(), 0);

    std::vector<toml::Table> fflocal, ffglobal, ffexternal;
    if(ff.count("local") == 1)
    {
        MJOLNIR_LOG_INFO("LocalForceField found");
        fflocal = toml::get<std::vector<toml::Table>>(
                toml_value_at(ff, "local", "[forcefields]"));
    }
    if(ff.count("global") == 1)
    {
        MJOLNIR_LOG_INFO("GlobalForceField found");
        ffglobal = toml::get<std::vector<toml::Table>>(
                toml_value_at(ff, "global", "[forcefields]"));
    }
    if(ff.count("external") == 1)
    {
        MJOLNIR_LOG_INFO("ExternalForceField found");
        ffexternal = toml::get<std::vector<toml::Table>>(
                toml_value_at(ff, "external", "[forcefields]"));
    }

    for(const auto kv: ff)
    {
        if(kv.first != "local" && kv.first != "global" && kv.first != "external")
        {
            std::cerr << "WARNING: unknown key `" << kv.first << "` appeared. ";
            std::cerr << "in [[forcefields]] table. It will be ignored.\n";
            MJOLNIR_LOG_WARN("unknown key ", kv.first, "appeared.");
        }
    }

    return ForceField<traitsT>(
            read_local_forcefield<traitsT>(std::move(fflocal)),
            read_global_forcefield<traitsT>(std::move(ffglobal)),
            read_external_forcefield<traitsT>(std::move(ffexternal)));
}


template<typename traitsT>
ForceField<traitsT>
read_forcefield(const toml::Table& data, std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_forcefield(), 0);

    const auto& ffs = toml_value_at(data, "forcefields", "<root>"
            ).cast<toml::value_t::Array>();
    if(ffs.size() <= N)
    {
        throw std::out_of_range("no enough forcefields: " + std::to_string(N));
    }
    MJOLNIR_LOG_INFO(ffs.size(), " forcefields are provided");
    MJOLNIR_LOG_INFO("using ", N, "-th forcefield");

    const auto& ff = ffs.at(N).cast<toml::value_t::Table>();
    if(ff.count("file_name") == 1)
    {
        MJOLNIR_SCOPE(ff.count("file_name") == 1, 1);
        if(ff.size() != 1)
        {
            std::cerr << "WARNING: [[forcefields]] has `file_name` key.";
            std::cerr << "When `file_name` is provided, all settings will be ";
            std::cerr << "read from the file, so other settings are ignored.\n";
            MJOLNIR_LOG_WARN("[[forcefields]] has file_name and other settings");
        }

        const std::string file_name = toml::get<std::string>(ff.at("file_name"));
        MJOLNIR_LOG_INFO("file_name = ", file_name);

        const auto forcefield_file = toml::parse(file_name);
        if(forcefield_file.count("forcefields") == 1)
        {
            MJOLNIR_LOG_INFO("`forcefields` value found in ", file_name);

            const auto forcefield_toml_type =
                forcefield_file.at("forcefields").type();
            if(forcefield_toml_type != toml::value_t::Table)
            {
                std::cerr << "FATAL: each [forcefields] should be provided as ";
                std::cerr << "a table in each file (" << file_name <<  ").\n";
                std::cerr << "     : because it represents one set of ";
                std::cerr << "forcefield. currently it is provided as ";
                std::cerr << forcefield_toml_type << ".\n";
                std::exit(1);
            }
            std::cerr << "WARNING: [forcefields] isn't necessary for splitted ";
            std::cerr << "toml file. you can define just [[local]], [[global]]";
            std::cerr << " and [[external]] forcefields in " << file_name << '\n';

            MJOLNIR_LOG_INFO("reading `forcefields` table");
            return read_forcefield_from_table<traitsT>(forcefield_file.at(
                    "forcefields").template cast<toml::value_t::Table>());
        }
        return read_forcefield_from_table<traitsT>(forcefield_file);
    }
    // else
    return read_forcefield_from_table<traitsT>(ff);
}


} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD
