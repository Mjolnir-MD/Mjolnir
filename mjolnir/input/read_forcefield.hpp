#ifndef MJOLNIR_READ_FORCEFIELD
#define MJOLNIR_READ_FORCEFIELD
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_interaction.hpp>

namespace mjolnir
{

template<typename traitsT>
LocalForceField<traitsT>
read_local_forcefield(std::vector<toml::Table> interactions,
                      const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_local_forcefield(), 0);
    MJOLNIR_LOG_INFO(interactions.size(),
                     " kinds of local interactions are found.");

    LocalForceField<traitsT> lff;
    for(const auto& interaction : interactions)
    {
        if(interaction.count("file_name") == 1)
        {
            MJOLNIR_SCOPE(interaction.count("file_name") == 1, 1);

            const std::string file_name = get_toml_value<std::string>(
                    interaction, "file_name", "[[local]]");
            MJOLNIR_LOG_INFO("file_name = ", file_name);

            if(interaction.size() != 1)
            {
                MJOLNIR_LOG_WARN(
                    "[[forcefields.local]] has `file_name` and other keys.");
                MJOLNIR_LOG_WARN(
                    "When `file_name` is provided, other values are ignored "
                    "because those are read from the specified file (",
                    file_name, ").");
            }

            const auto forcefield_file = toml::parse(input_path + file_name);
            if(forcefield_file.count("forcefields") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[local] should be provided as a root object of file ",
                    file_name, ". but [[forcefields]] table found");
                std::exit(1);
            }
            if(forcefield_file.count("local") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[local] should be provided as a root object of file ",
                    file_name, ". but [local] table found");

                if(forcefield_file.at("local").type() != toml::value_t::Table)
                {
                    MJOLNIR_LOG_ERROR("type of `local` is different from toml::"
                                      "Table in file (", file_name, ").");
                    MJOLNIR_LOG_ERROR("note: [[...]] means Array-of-Tables. "
                                      "please take care.");
                    std::exit(1);
                }
                lff.emplace(read_local_interaction<traitsT>(
                    get_toml_value<toml::Table>(
                        forcefield_file, "local", file_name)));
            }
            else
            {
                lff.emplace(read_local_interaction<traitsT>(forcefield_file));
            }
        }
        else
        {
            lff.emplace(read_local_interaction<traitsT>(interaction));
        }
    }
    return lff;
}

template<typename traitsT>
GlobalForceField<traitsT>
read_global_forcefield(std::vector<toml::Table> interactions,
                       const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_forcefield(), 0);
    MJOLNIR_LOG_INFO(interactions.size(),
                     " kinds of local interactions are found.");

    GlobalForceField<traitsT> gff;
    for(const auto& interaction : interactions)
    {
        if(interaction.count("file_name") == 1)
        {
            MJOLNIR_SCOPE(interaction.count("file_name") == 1, 1);

            const std::string file_name = get_toml_value<std::string>(
                    interaction, "file_name", "[[global]]");
            MJOLNIR_LOG_INFO("file_name = ", file_name);

            if(interaction.size() != 1)
            {
                MJOLNIR_LOG_WARN(
                    "[[forcefields.global]] has `file_name` and other keys.");
                MJOLNIR_LOG_WARN(
                    "When `file_name` is provided, other values are ignored "
                    "because those are read from the specified file (",
                    file_name, ").");
            }

            const auto forcefield_file = toml::parse(input_path + file_name);
            if(forcefield_file.count("forcefields") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[global] should be provided as a root object of file ",
                    file_name, ". but [[forcefields]] table found");
                std::exit(1);
            }
            if(forcefield_file.count("global") == 1)
            {
                MJOLNIR_LOG_ERROR(
                    "[global] should be provided as a root object of file ",
                    file_name, ". but [global] key found");

                if(forcefield_file.at("global").type() != toml::value_t::Table)
                {
                    MJOLNIR_LOG_ERROR("type of `global` is different from "
                                      "toml::Table in file (", file_name, ").");
                    MJOLNIR_LOG_ERROR("note: [[...]] means Array-of-Tables. "
                                      "please take care.");
                    std::exit(1);
                }
                gff.emplace(read_global_interaction<traitsT>(
                    get_toml_value<toml::Table>(
                        forcefield_file, "global", file_name)));
            }
            else
            {
                gff.emplace(read_global_interaction<traitsT>(forcefield_file));
            }
        }
        else
        {
            gff.emplace(read_global_interaction<traitsT>(interaction));
        }
    }
    return gff;
}

template<typename traitsT>
ExternalForceField<traitsT>
read_external_forcefield(std::vector<toml::Table> interactions,
                         const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_external_forcefield(), 0);
    MJOLNIR_LOG_INFO(interactions.size(),
                     " kinds of external interactions are found.");

    ExternalForceField<traitsT> eff;
    for(const auto& interaction : interactions)
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

    const auto& files = get_toml_value<toml::Table>(data, "files", "<root>");
    std::string input_path_("./");
    if(files.count("input_path") == 1)
    {
        input_path_ = get_toml_value<std::string>(files, "input_path", "[files]");
    }
    const auto input_path(input_path_); // add constness
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
#endif// MJOLNIR_READ_FORCEFIELD
