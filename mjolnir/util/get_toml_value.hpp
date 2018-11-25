#ifndef MJOLNIR_UTIL_GET_TOML_VALUE
#define MJOLNIR_UTIL_GET_TOML_VALUE
#include <extlib/toml/toml/toml.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <string>
#include <sstream>
#include <stdexcept>

namespace mjolnir
{

// utility functions to get a value from toml::table.
// To show the informative error message, use it instead of the raw toml::get.

template<typename T>
decltype(toml::get<T>(std::declval<toml::value>()))
get_toml_value(const toml::table& tab, const std::string& key,
               const std::string& tablename)
{
    try
    {
        return toml::get<T>(tab.at(key));
    }
    catch(const std::out_of_range& oor)
    {
        throw_exception<std::runtime_error>(
            "mjolnir: while reading toml file: key(", key,
            ") in table(", tablename, ") missing");
    }
    catch(const toml::bad_get& bg)
    {
        throw_exception<std::runtime_error>(
            "mjolnir: while reading toml file: key(", key,
            ") in table(", tablename, "): ", bg.what());
    }
}

template<typename T>
decltype(toml::get<T>(std::declval<toml::value>()))
get_toml_value(const toml::table& tab, std::initializer_list<std::string> keys,
               const std::string& tablename)
{
    for(const auto& key : keys)
    {
        const auto iter = tab.find(key);
        if(iter != tab.end())
        {
            try
            {
                return toml::get<T>(iter->second);
            }
            catch(const toml::bad_get& bg)
            {
                throw_exception<std::runtime_error>(
                    "mjolnir: while reading toml file: key(", key,
                    ") in table(", tablename, "): ", bg.what());
            }
        }
    }

    // not found. throw an error.

    std::ostringstream oss;
    oss << "mjolnir: while reading toml file: none of the key(";
    for(const auto& key : keys) {oss << key << ", ";}
    oss << ") is found in table (" << tablename << ").";
    throw std::runtime_error(oss.str());
}

} // mjolnir
#endif// MJOLNIR_INPUT_GET_TOML_VALUE
