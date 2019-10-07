#ifndef MJOLNIR_INPUT_READ_TABLE_FROM_FILE_HPP
#define MJOLNIR_INPUT_READ_TABLE_FROM_FILE_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/format_nth.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/input/utility.hpp>

// this function is used by read_system and read_forcefield.

namespace mjolnir
{
// ---------------------------------------------------------------------------
// a function to read just one table (normal case).

template<typename C, template<typename...> class T, template<typename...> class A>
toml::basic_value<C, T, A>
read_table_from_file(const toml::basic_value<C, T, A>& root,
        const std::string& name, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_NOTICE("reading [[", name, "]].");

    const auto& tables = toml::find(root, name).as_array();
    if(tables.empty())
    {
        throw_exception<std::runtime_error>(
            "[error] mjolnir::read_table_from_file: no ", name, " defined.");
    }

    const auto& table = tables.front();
    if(table.as_table().count("file_name") == 0)
    {
        // no "file_name" is defined. we don't need any file.
        return table;
    }

    // file_name is provided. we need to read it.
    const auto file_name = input_path + toml::find<std::string>(table, "file_name");

    MJOLNIR_LOG_NOTICE("[[", name, "]] is defined in ", file_name);
    if(table.as_table().size() != 1)
    {
        MJOLNIR_LOG_WARN("When `file_name` is provided, all other keys will be"
                         "ignored.");
        check_keys_available(table, {"file_name"_s});
    }

    MJOLNIR_LOG_NOTICE("reading ", file_name, " ...");
    const auto table_file = toml::parse(file_name);
    MJOLNIR_LOG_NOTICE(" done.");

    if(table_file.as_table().count(name) == 0)
    {
        throw_exception<std::out_of_range>("[error] mjolnir::"
            "read_table_from_file: table [[", name, "]] not found in toml file"
            "\n --> ", file_name, "\n | the file should define [[", name,
            "]] table and define values in it.");
    }

    if(table_file.as_table().at(name).is_array())
    {
        return toml::find(table_file, name).as_array().front();
    }
    return toml::find(table_file, name);
}

// ---------------------------------------------------------------------------
// a function to read N-th table (multiple forcefield/system case).


template<typename C, template<typename...> class T, template<typename...> class A>
toml::basic_value<C, T, A>
read_table_from_file(const toml::basic_value<C, T, A>& root,
    const std::string& name, const std::size_t N, const std::string& input_path)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_NOTICE("reading ", format_nth(N), " [[", name, "]].");

    const auto& tables = toml::find(root, name).as_array();
    if(tables.size() <= N)
    {
        throw_exception<std::runtime_error>("[error] mjolnir::"
            "read_table_from_file: no enough number of ", name, " provided.");
    }
    MJOLNIR_LOG_INFO(tables.size(), " tables are provided");
    MJOLNIR_LOG_INFO("using ", N, "-th table");

    const auto& table = tables.at(N);
    if(table.as_table().count("file_name") == 0)
    {
        // no "file_name" is defined. we don't need any file.
        return table;
    }

    // file_name is provided. we need to read it.
    const auto file_name = input_path + toml::find<std::string>(table, "file_name");

    MJOLNIR_LOG_NOTICE("[[", name, "]] is defined in ", file_name);
    if(table.as_table().size() != 1)
    {
        MJOLNIR_LOG_WARN("When `file_name` is provided, all other keys will be"
                         "ignored.");
        check_keys_available(table, {"file_name"_s});
    }

    MJOLNIR_LOG_NOTICE("reading ", file_name, " ...");
    const auto table_file = toml::parse(file_name);
    MJOLNIR_LOG_NOTICE(" done.");

    if(table_file.as_table().count(name) == 0)
    {
        throw_exception<std::out_of_range>("[error] mjolnir::"
            "read_table_from_file: table [[", name, "]] not found in toml file"
            "\n --> ", file_name, "\n | the file should define [[", name,
            "]] table and define values in it.");
    }

    if(table_file.as_table().at(name).is_array())
    {
        return toml::find(table_file, name).as_array().front();
    }
    return toml::find(table_file, name);
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>
read_table_from_file(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>& root,
                     const std::string& name, const std::string& input_path);
extern template toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>
read_table_from_file(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>& root,
                     const std::string& name, const std::size_t N, const std::string& input_path);
#endif

} // mjolnir
#endif//MJOLNIR_INPUT_READ_TABLE_FROM_FILE_HPP
