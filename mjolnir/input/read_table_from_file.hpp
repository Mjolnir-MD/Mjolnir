#ifndef MJOLNIR_INPUT_READ_TABLE_FROM_FILE_HPP
#define MJOLNIR_INPUT_READ_TABLE_FROM_FILE_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/format_nth.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/input/utility.hpp>

// some tables such as [simulator], [[systems]], [[forcefields]], and
// [[forcefields.local]] are allowed to be provided as another file in the
// following way.
//
// ```toml:main.toml
// [[systems]]
// file_name = "initial_structure.toml"
// ```
//
// ```toml:initial_structure.toml
// [[systems]]
// attributes.temperature = 300.0
// particles = [
// # ...
// ]
// ```
//
// To allow both, this function checks if the table has `file_name` and parse
// the file if the file_name is provided. The path is under `files.input.path`.

namespace mjolnir
{

template<typename C, template<typename...> class T, template<typename...> class A>
toml::basic_value<C, T, A>
read_table_from_file(const toml::basic_value<C, T, A>& target,
                     const std::string& table_name, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if(target.as_table().count("file_name") == 0)
    {
        // no file_name = "..." is defined. we don't need to parse another file
        return target;
    }
    if(target.as_table().size() != 1)
    {
        MJOLNIR_LOG_WARN("When `file_name` is provided, all other keys will be"
                         "ignored.");
        check_keys_available(target, {"file_name"_s});
    }

    // file_name is provided. we need to read it.
    const auto file_name = input_path + toml::find<std::string>(target, "file_name");
    MJOLNIR_LOG_NOTICE("table ", table_name, " is defined in ", file_name);

    MJOLNIR_LOG_NOTICE("reading ", file_name, " ...");
    const auto table_from_file = toml::parse(file_name);
    MJOLNIR_LOG_NOTICE("done.");

    if(table_from_file.as_table().count(table_name) == 0)
    {
        throw_exception<std::out_of_range>("[error] mjolnir::"
            "read_table_from_file: table [", table_name, "] not found in toml"
            " file\n --> ", file_name, "\n  | the file should define [",
            table_name, "] table and define values in it.");
    }
    // if it defines an array of tables, use the first one.
    //
    // ```toml:initial_structure.toml
    // [[systems]] # this is an array of table
    // particles = [
    // # ...
    // ]
    if(table_from_file.as_table().at(table_name).is_array())
    {
        return toml::find(table_from_file, table_name).as_array().front();
    }
    return toml::find(table_from_file, table_name);
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>
read_table_from_file(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>& root,
                     const std::string& table_name, const std::string& input_path);
#endif

} // mjolnir
#endif//MJOLNIR_INPUT_READ_TABLE_FROM_FILE_HPP
