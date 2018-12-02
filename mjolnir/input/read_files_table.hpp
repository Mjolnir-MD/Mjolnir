#ifndef MJOLNIR_READ_FILES_TABLE_H
#define MJOLNIR_READ_FILES_TABLE_H
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// this function may be callen from other read_* functions.
inline std::string read_input_path(const toml::Table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_input_path(), 0);

    const auto& files = get_toml_value<toml::Table>(root, "files", "<root>");

    std::string input_path;
    if(files.count("input") == 1)
    {
        const auto& input = get_toml_value<toml::Table>(files, "input", "[files]");
        if(input.count("path") == 1)
        {
            input_path = get_toml_value<std::string>(input, "path", "[files.input]");
        }
    }
    return input_path;
}

} // mjolnir
#endif// MJOLNIR_READ_FILES_TABLE_H
