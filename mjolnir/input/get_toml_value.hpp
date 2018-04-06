#ifndef MJOLNIR_INPUT_GET_TOML_VALUE
#define MJOLNIR_INPUT_GET_TOML_VALUE
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <string>
#include <stdexcept>

namespace mjolnir
{

inline toml::value const&
toml_value_at(const toml::Table& tab, const std::string& key,
              const std::string& tablename)
{
    try
    {
        return tab.at(key);
    }
    catch(const std::out_of_range& oor)
    {
        throw_exception<std::runtime_error>(
            "mjolnir: while reading toml file: key(", key,
            ") in table(", tablename, ") missing");
    }
}

} // mjolnir
#endif// MJOLNIR_INPUT_GET_TOML_VALUE
