#ifndef JARNGREIPR_WRITE_TOML_VALUE_HPP
#define JARNGREIPR_WRITE_TOML_VALUE_HPP
#include <extlib/toml/toml.hpp>

namespace jarngreipr
{

template<typename charT, typename traits>
void write_toml_value(
        std::basic_ostream<charT, traits>& os, const toml::value& v)
{
    switch(v.type())
    {
        case toml::value_t::Boolean:
        {
            os << std::boolalpha << toml::get<bool>(v);
            return;
        }
        case toml::value_t::Integer:
        {
            os << toml::get<toml::Integer>(v);
            return;
        }
        case toml::value_t::Floating:
        {
            os << std::fixed << toml::get<toml::Floating>(v);
            return;
        }
        case toml::value_t::String:
        {
            os << '"' << toml::get<toml::Floating>(v) << '"';
            return;
        }
        case toml::value_t::Array:
        {
            os << '[';
            for(const auto& item : v.cast<toml::value_t::Array>())
            {
                write_toml_value(os, item);
                os << ',';
            }
            os << ']';
            return;
        }
        case toml::value_t::Table:
        {
            os << '{';
            for(const auto& kv : v.cast<toml::value_t::Table>())
            {
                os << kv.first << '=';
                write_toml_value(os, kv.second);
                os << ',';
            }
            os << '}';
            return;
        }
        default:
            throw std::runtime_error("write_toml_value: not a simple value: " +
                    toml::stringize(v.type()));
    }
}

} // jarngreipr
#endif// JARNGREIPR_WRITE_TOML_VALUE_HPP
