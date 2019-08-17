#ifndef MJOLNIR_INPUT_READ_PATH_HPP
#define MJOLNIR_INPUT_READ_PATH_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// this function may be callen from other read_* functions.
inline std::string read_input_path(const toml::value& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto& files = toml::find(root, "files");

    std::string input_path("./");
    if(files.as_table().count("input") == 1)
    {
        const auto& input = toml::find(files, "input");
        if(input.as_table().count("path") == 1)
        {
            input_path = toml::find<std::string>(input, "path");
            if(input_path.back() != '/') {input_path += '/';}
        }
    }
    return input_path;
}

inline std::string read_output_path(const toml::value& root)
{
    // This function does not output log because it might be called when
    // the logger has not been initialized yet.
    // Logging should be done after calling this function.

    const auto& files = toml::find(root, "files");

    std::string output_path("./");
    if(files.as_table().count("output") == 1)
    {
        const auto& output = toml::find(files, "output");
        if(output.as_table().count("path") == 1)
        {
            output_path = toml::find<std::string>(output, "path");
            if(output_path.back() != '/') {output_path += '/';}
        }
    }
    return output_path;
}

} // mjolnir
#endif// MJOLNIR_READ_PATH_HPP
