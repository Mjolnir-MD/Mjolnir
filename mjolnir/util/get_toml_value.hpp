#ifndef MJOLNIR_UTIL_GET_TOML_VALUE
#define MJOLNIR_UTIL_GET_TOML_VALUE
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/type_traits.hpp>
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

// one of the toml values. needs no conversion. so it can return a reference.
template<typename T, typename std::enable_if<disjunction<
    std::is_same<T, toml::Boolean >,
    std::is_same<T, toml::Integer >,
    std::is_same<T, toml::Float   >,
    std::is_same<T, toml::String  >,
    std::is_same<T, toml::Datetime>,
    std::is_same<T, toml::Array   >,
    std::is_same<T, toml::Table   >
    >::value, std::nullptr_t>::type = nullptr>
T const& get_toml_value(const toml::Table& tab, const std::string& key,
                        const std::string& tablename)
{
    try
    {
        return tab.at(key).cast<toml::value_traits<T>::type_index>();
    }
    catch(const std::out_of_range& oor)
    {
        throw_exception<std::runtime_error>(
            "mjolnir: while reading toml file: key(", key,
            ") in table(", tablename, ") missing");
    }
    catch(const toml::type_error& tte)
    {
        throw_exception<std::runtime_error>(
            "mjolnir: while reading toml file: key(", key,
            ") in table(", tablename, ") has type `",
            tab.at(key).type(), "`, expected `",
            toml::value_traits<T>::type_index, "`.");
    }
}

// none of the toml values. needs conversion.
template<typename T, typename std::enable_if<conjunction<
    negation<std::is_same<T, toml::Boolean >>,
    negation<std::is_same<T, toml::Integer >>,
    negation<std::is_same<T, toml::Float   >>,
    negation<std::is_same<T, toml::String  >>,
    negation<std::is_same<T, toml::Datetime>>,
    negation<std::is_same<T, toml::Array   >>,
    negation<std::is_same<T, toml::Table   >>
    >::value, std::nullptr_t>::type = nullptr>
T get_toml_value(const toml::Table& tab, const std::string& key,
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
    catch(const toml::type_error& tte)
    {
        throw_exception<std::runtime_error>(
            "mjolnir: while reading toml file: key(", key,
            ") in table(", tablename, ") has type `",
            tab.at(key).type(), "`, expected `",
            toml::value_traits<T>::type_index, "`.");
    }
}

// one of the toml values. needs no conversion. so it can return a reference.
template<typename T, typename std::enable_if<disjunction<
    std::is_same<T, toml::Boolean >,
    std::is_same<T, toml::Integer >,
    std::is_same<T, toml::Float   >,
    std::is_same<T, toml::String  >,
    std::is_same<T, toml::Datetime>,
    std::is_same<T, toml::Array   >,
    std::is_same<T, toml::Table   >
    >::value, std::nullptr_t>::type = nullptr>
T const& get_toml_value(const toml::Table& tab,
                        std::initializer_list<std::string> keys,
                        const std::string& tablename)
{
    for(const auto& key : keys)
    {
        const auto iter = tab.find(key);
        if(iter != tab.end())
        {
            return tab.at(key).cast<toml::value_traits<T>::type_index>();
        }
    }

    std::ostringstream oss;
    oss << "mjolnir: while reading toml file: none of the key(";
    for(const auto& key : keys)
    {
        oss << key << ", ";
    }
    oss << ") are found in table (" << tablename << ").";
    throw std::runtime_error(oss.str());
}

// none of the toml values. needs conversion.
template<typename T, typename std::enable_if<conjunction<
    negation<std::is_same<T, toml::Boolean >>,
    negation<std::is_same<T, toml::Integer >>,
    negation<std::is_same<T, toml::Float   >>,
    negation<std::is_same<T, toml::String  >>,
    negation<std::is_same<T, toml::Datetime>>,
    negation<std::is_same<T, toml::Array   >>,
    negation<std::is_same<T, toml::Table   >>
    >::value, std::nullptr_t>::type = nullptr>
T get_toml_value(const toml::Table& tab,
                 std::initializer_list<std::string> keys,
                 const std::string& tablename)
{
    for(const auto& key : keys)
    {
        const auto iter = tab.find(key);
        if(iter != tab.end())
        {
            return toml::get<T>(tab.at(key));
        }
    }

    std::ostringstream oss;
    oss << "mjolnir: while reading toml file: none of the key(";
    for(const auto& key : keys)
    {
        oss << key << ", ";
    }
    oss << ") are found in table (" << tablename << ").";
    throw std::runtime_error(oss.str());}


} // mjolnir
#endif// MJOLNIR_INPUT_GET_TOML_VALUE
