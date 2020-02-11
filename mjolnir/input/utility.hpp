#ifndef MJOLNIR_INPUT_UTILITY_HPP
#define MJOLNIR_INPUT_UTILITY_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/type_traits.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/math/Matrix.hpp>

namespace mjolnir
{

// This check all the keys in a table are found in a list.
//     If there is a key that is not found in the range, it warns about the
// corresponding value will be ignored.
//     In order to allow optional keys, it only checks all the keys are found
// in the specified container.
//
// Use it as the following.
// ```cpp
// check_keys_available(table, {"foo"_s, "bar"_s, "baz"_s});
// ```
inline bool check_keys_available(const toml::value& table,
                                 std::initializer_list<std::string> list)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    // no logger scope here. use the parent scope.

    bool all_available = true;
    for(const auto& kv : table.as_table())
    {
        if(list.end() == std::find(list.begin(), list.end(), kv.first))
        {
            std::ostringstream oss;
            oss << "unknown value \"" << kv.first << "\" found. this "
                << kv.second.type() << " will never be used.";
            const auto err_msg = toml::format_error(
                    oss.str(), kv.second, "this will be ignored");
            // workaround to skip auto-added [error].
            // Normally this function is called only when the input file
            // contains an invalid key. So it should not be a hotspot.
            MJOLNIR_LOG_WARN(err_msg.substr(err_msg.find(oss.str())));
            all_available = false;
        }
    }
    return all_available;
}

template<typename T>
typename std::enable_if<negation<std::is_same<T, std::string>>::value, T>::type
find_parameter(const toml::value& params, const toml::value& env,
               const std::string& name)
{
    using ::mjolnir::literals::string_literals::operator"" _s;
    if(!params.is_table() || params.as_table().count(name) != 1)
    {
        throw std::out_of_range(toml::format_error("[error] value "_s + name +
            " does not exists"_s, params, "in this table"_s));
    }
    const toml::value& p = params.as_table().at(name);
    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_uninitialized())
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" used but no env is defined"_s, params, "used here"_s));
        }
        if(!env.is_table() || env.as_table().count(var) != 1)
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" does not exists"_s, env, "in this table"_s));
        }
        return toml::get<T>(env.as_table().at(var));
    }
    return toml::get<T>(p);
}

// If the expected value is std::string, it is difficult to distinguish
// variable name and the value itself. In that case, env would just be ignored.
template<typename T>
typename std::enable_if<std::is_same<T, std::string>::value, std::string>::type
find_parameter(const toml::value& params, const toml::value& /* env */,
               const std::string& name)
{
    if(!params.is_table() || params.as_table().count(name) != 1)
    {
        throw std::out_of_range(toml::format_error("[error] value \""_s + name +
            "\" does not exists"_s, params, "in this table"_s));
    }
    const toml::value& p = params.as_table().at(name);
    return toml::get<T>(p);
}

// The current version of Mjolnir allow to use unicode character when defining
// a parameter. But since it is a bit ambiguous because of some reasons. E.g.
// - The same character might appear in the unicode table several time with
//   different styles
// - When the different parameters appear in different names, it is ofcourse
//   ambiguous which one to be used.
//
// So I decided to deprecate them. At first I thought that a unicode names
// ware covenient to reduce the size of input file. But now, since env is
// introduced, it is more effective for reducing the file size. So if it works
// nice after the next release, I will deprecate unicode names and warn if used.
//
// This function was introduced to support those
// multi-named parameters but is planned to be removed in the later release.
template<typename T>
typename std::enable_if<negation<std::is_same<T, std::string>>::value, T>::type
find_parameter(const toml::value& params, const toml::value& env,
               const std::string& name1,  const std::string& name2)
{
    if(!params.is_table() || (params.as_table().count(name1) == 0 &&
                              params.as_table().count(name2) == 0))
    {
        throw std::out_of_range(toml::format_error("[error] value \""_s + name1 +
            "\" or \""_s + name2 + "\" does not exists"_s, params, "in this table"_s));
    }
    // name1 has priority.
    const toml::value& p = (params.as_table().count(name1) == 1) ?
                            params.as_table().at(name1)          :
                            params.as_table().at(name2)          ;
    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_uninitialized())
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" used but no env is defined"_s, params, "used here"_s));
        }
        if(!env.is_table() || env.as_table().count(var) != 1)
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" does not exists"_s, env, "in this table"_s));
        }
        return toml::get<T>(env.as_table().at(var));
    }
    return toml::get<T>(p);
}

} // mjolnir

namespace toml
{

// enable to get mjolnir::Vector as the following.
// ```cpp
// const auto v = toml::find<mjolnir::Vector<double, 3>>(table, "position");
// ```

template<>
struct from<mjolnir::math::Matrix<double, 3, 1>>
{
    template<typename C, template<typename ...> class M,
             template<typename ...> class A>
    static mjolnir::math::Matrix<double, 3, 1> from_toml(const basic_value<C, M, A>& v)
    {
        return mjolnir::math::Matrix<double, 3, 1>(get<std::array<double, 3>>(v));
    }
};

template<>
struct from<mjolnir::math::Matrix<float, 3, 1>>
{
    template<typename C, template<typename ...> class M,
             template<typename ...> class A>
    static mjolnir::math::Matrix<float, 3, 1> from_toml(const basic_value<C, M, A>& v)
    {
        return mjolnir::math::Matrix<float, 3, 1>(get<std::array<float, 3>>(v));
    }
};

} // toml
#endif// MJOLNIR_INPUT_UTILITY_HPP
