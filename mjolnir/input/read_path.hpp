#ifndef MJOLNIR_INPUT_READ_FILES_TABLE_HPP
#define MJOLNIR_INPUT_READ_FILES_TABLE_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// this function may be callen from other read_* functions.
inline std::string read_input_path(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto& files = toml::find<toml::table>(root, "files");

    std::string input_path("./");
    if(files.count("input") == 1)
    {
        const auto& input = toml::find<toml::table>(files, "input");
        if(input.count("path") == 1)
        {
            input_path = toml::find<std::string>(input, "path");
            if(input_path.back() != '/') {input_path += '/';}
        }
    }
    return input_path;
}

inline std::string read_output_path(const toml::table& root)
{
    // This function does not output log because it might be called when
    // the logger has not been initialized yet.
    // Logging should be done after calling this function.

    const auto& files = toml::find<toml::table>(root, "files");

    std::string output_path("./");
    if(files.count("output") == 1)
    {
        const auto& output = toml::find<toml::table>(files, "output");
        if(output.count("path") == 1)
        {
            output_path = toml::find<std::string>(output, "path");
            if(output_path.back() != '/') {output_path += '/';}
        }
    }
    return output_path;
}

} // mjolnir
#endif// MJOLNIR_READ_FILES_TABLE_H
