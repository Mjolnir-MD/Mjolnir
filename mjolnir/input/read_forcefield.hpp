#ifndef MJOLNIR_READ_FORCEFIELD_HPP
#define MJOLNIR_READ_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_forcefield.hpp>
#include <mjolnir/input/read_global_forcefield.hpp>
#include <mjolnir/input/read_external_forcefield.hpp>
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
ForceField<traitsT>
read_forcefield_from_table(const toml::Table& ff, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_forcefield_from_table(), 0);

    std::vector<toml::Table> fflocal, ffglobal, ffexternal;
    if(ff.count("local") == 1)
    {
        MJOLNIR_LOG_INFO("LocalForceField found");
        fflocal = get_toml_value<std::vector<toml::Table>>(
                ff, "local", "[[forcefields]]");
    }
    if(ff.count("global") == 1)
    {
        MJOLNIR_LOG_INFO("GlobalForceField found");
        ffglobal = get_toml_value<std::vector<toml::Table>>(
                ff, "global", "[[forcefields]]");
    }
    if(ff.count("external") == 1)
    {
        MJOLNIR_LOG_INFO("ExternalForceField found");
        ffexternal = get_toml_value<std::vector<toml::Table>>(
                ff, "external", "[[forcefields]]");
    }

    for(const auto kv: ff)
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
read_forcefield(const toml::Table& data, std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_forcefield(), 0);
    MJOLNIR_LOG_NOTICE("reading ", N, "-th [[forcefield]].");

    const auto input_path = read_input_path(data);
    MJOLNIR_LOG_INFO("input path is -> ", input_path);

    const auto& ffs = get_toml_value<toml::Array>(data, "forcefields", "<root>");
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

        const std::string file_name =
            get_toml_value<std::string>(ff, "file_name", "[[forcefields]]");
        MJOLNIR_LOG_INFO("file_name = ", file_name);

        if(ff.size() != 1)
        {
            MJOLNIR_LOG_WARN("[[forcefields]] has `file_name` and other keys.");
            MJOLNIR_LOG_WARN("When `file_name` is provided, other values are "
                             "ignored because those are read from the specified"
                             " file (", file_name, ").");
        }

        MJOLNIR_LOG_NOTICE("forcefield is defined in `", input_path, file_name, "`.");
        const auto forcefield_file = toml::parse(input_path + file_name);
        if(forcefield_file.count("forcefields") == 1)
        {
            MJOLNIR_LOG_WARN("in `forcefield` file, root object is treated as "
                             "one of the [[forcefields]] tables.");
            MJOLNIR_LOG_WARN("define just [[local]], [[global]] and "
                             "[[external]] forcefields in ", file_name);
            MJOLNIR_LOG_WARN("in ", file_name, ", [forcefields] key found."
                             "trying to read it as a forcefield setup.");

            if(forcefield_file.at("forcefields").type() != toml::value_t::Table)
            {
                MJOLNIR_LOG_ERROR("type of `forcefields` is different from "
                                  "toml::Table in file (", file_name, ").");
                MJOLNIR_LOG_ERROR("note: [[...]] means Array-of-Tables. "
                                  "please take care.");
                std::exit(1);
            }
            return read_forcefield_from_table<traitsT>(forcefield_file.at(
                    "forcefields").template cast<toml::value_t::Table>(),
                    input_path);
        }
        return read_forcefield_from_table<traitsT>(forcefield_file, input_path);
    }

    // all-in-one file
    return read_forcefield_from_table<traitsT>(ff, input_path);
}

} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD_HPP
