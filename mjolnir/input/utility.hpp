#ifndef MJOLNIR_INPUT_UTILITY_HPP
#define MJOLNIR_INPUT_UTILITY_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/util/type_traits.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/math/Matrix.hpp>

namespace mjolnir
{

// find_parameter is a utility function to support the following functionality.
//
// ```toml
// [[forcefields.global]]
// env.sigma = 1.0
// parameters = [
//     {index = 1, sigma = "sigma"},
//     {index = 2, sigma = 2.0},
//     # ...
// ]
// ```
//
// First, it searches `params` with the `name`. If the corresponding value is a
// string, it searches `env` with the string.
template<typename T>
T find_parameter(const toml::value& params, const toml::value& env,
                 const std::string& name)
{
    static_assert(!std::is_same<T, std::string>::value,
                  "string value cannot be aliased");

    using ::mjolnir::literals::string_literals::operator"" _s;
    if(!params.is_table() || !params.contains(name))
    {
        throw std::out_of_range(toml::format_error("[error] value "_s + name +
            " does not exists"_s, params, "in this table"_s));
    }
    const toml::value& p = params.at(name);
    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_uninitialized())
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" used but no env is defined"_s, params, "used here"_s));
        }
        if(!env.is_table() || !env.contains(var))
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" does not exists"_s, env, "in this table"_s));
        }
        return toml::find<T>(env, var);
    }
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
    if(!params.is_table() || (!params.contains(name1) && !params.contains(name2)))
    {
        throw std::out_of_range(toml::format_error("[error] value \""_s + name1 +
            "\" or \""_s + name2 + "\" does not exists"_s, params, "in this table"_s));
    }
    // name1 has priority.
    const toml::value& p = params.contains(name1) ? params.at(name1) : params.at(name2);
    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_uninitialized())
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" used but no env is defined"_s, params, "used here"_s));
        }
        if(!env.is_table() || !env.contains(var))
        {
            throw std::out_of_range(toml::format_error("[error] named variable \""_s +
                var + "\" does not exists"_s, env, "in this table"_s));
        }
        return toml::find<T>(env, var);
    }
    return toml::get<T>(p);
}

// find_parameter with an optional value, opt.
template<typename T>
T find_parameter_or(const toml::value& params, const toml::value& env,
                    const std::string& name, const T& opt) noexcept
{
    static_assert(!std::is_same<T, std::string>::value,
                  "string value cannot be aliased");

    using ::mjolnir::literals::string_literals::operator"" _s;
    if(!params.is_table() || !params.contains(name))
    {
        return opt;
    }
    const toml::value& p = params.at(name);
    if(p.is_string())
    {
        return toml::find_or(env, p.as_string(), opt);
    }
    return toml::get_or(p, opt);
}

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
        if(kv.first == "include") {continue;}
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

// It checks that the index in the parameters does not appear more than once.
//     If the index appear twice, it raise the error represent the line of
// input file define the index.
// ```toml
// parameters = [
//     {index =  0, offset = 10, radius = 2.0},
//     {index = 10, radius = 3.0}, # <- overlap! parameter value is ambiguous.
// ]
// ```
template<typename parameterT>
void check_parameter_overlap(const toml::value& env, const toml::array& setting,
        std::vector<std::pair<std::size_t, parameterT>>& parameters)
{
    if(parameters.empty()) {return ;}

    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using value_type = std::pair<std::size_t, parameterT>;

    // sort the parameters by its index.
    std::sort(parameters.begin(), parameters.end(),
            [](const value_type& lhs, const value_type& rhs) noexcept -> bool {
                return lhs.first < rhs.first;
            });
    // Since the parameters are already sorted, all the overlapping parameters
    // should be next to each other. We can find it by adjacent_find.
    const auto overlap = std::adjacent_find(parameters.begin(), parameters.end(),
            [](const value_type& lhs, const value_type& rhs) noexcept -> bool {
                return lhs.first == rhs.first;
            });

    // If overlap is found, generate the error message and throw an error.
    if(overlap != parameters.end())
    {
        const std::size_t overlapped_idx = overlap->first;
        MJOLNIR_LOG_ERROR("parameter for ", overlapped_idx, " defined twice");

        // define a comparator to find the overlapping two parameters from
        // the toml array.
        const auto overlap_finder =
            [overlapped_idx, &env](const toml::value& v) -> bool {
                const auto ofs = find_parameter_or<std::int64_t>(v, env, "offset", 0);
                const auto idx = find_parameter   <std::size_t >(v, env, "index") + ofs;
                return idx == overlapped_idx;
            };

        const auto overlapped1 =
            std::find_if(setting.begin(), setting.end(), overlap_finder);
        assert(overlapped1 != setting.end());

        const auto overlapped2 =
            std::find_if(std::next(overlapped1), setting.end(), overlap_finder);
        assert(overlapped2 != setting.end());

        throw_exception<std::runtime_error>(toml::format_error(
            "duplicate parameter definitions",
            *overlapped1, "this defined twice", *overlapped2, "here"));
    }
    return ;
}

//
// if the file extension is the same as expected, return true.
// `expected` should contain the dot. e.g. expected = ".xyz"
//
inline bool file_extension_is(const std::string& filename,
                              const std::string& expected)
{
    if(filename.size() < expected.size())
    {
        return false;
    }
    const auto last_dot = filename.find_last_of('.');
    if(last_dot == std::string::npos)
    {
        return false;
    }
    return filename.substr(last_dot) == expected;
}

inline void merge_toml_tables(toml::value& table, const toml::value& other)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using ::mjolnir::literals::string_literals::operator"" _s;
    assert(table.is_table());
    assert(other.is_table());

    for(const auto& kv : other.as_table())
    {
        if(table.contains(kv.first))
        {
            if(table.at(kv.first).is_table())
            {
                MJOLNIR_LOG_INFO("merging table, ", kv.first);
                merge_toml_tables(table.at(kv.first), kv.second);
            }
            else
            {
                throw std::runtime_error(toml::format_error("value \""_s +
                    kv.first + "\" duplicates."_s, table, "first defined here"_s,
                    kv.second, "and also defined here"_s));
            }
        }
        else
        {
            MJOLNIR_LOG_INFO("inserting new value, ", kv.first);
            table.as_table().emplace(kv.first, kv.second);
        }
    }
    return;
}

inline void expand_include(toml::value& v)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    if(!v.is_table() && !v.is_array()) {return;}

    if(v.is_table())
    {
        for(auto& kv : v.as_table())
        {
            if(kv.first == "include") {continue;}
            expand_include(kv.second);
        }

        // expand include in this table
        if(v.contains("include"))
        {
            MJOLNIR_LOG_INFO("\"include\" found.");
            const auto& input_path = get_input_path();
            if(v.at("include").is_array())
            {
                for(auto fname : toml::find<std::vector<std::string>>(v, "include"))
                {
                    MJOLNIR_LOG_INFO("expanding file ", input_path, fname);
                    merge_toml_tables(v, toml::parse(input_path + fname));
                }
            }
            else
            {
                const auto& fname = toml::find<std::string>(v, "include");
                MJOLNIR_LOG_INFO("expanding file ", input_path, fname);
                merge_toml_tables(v, toml::parse(input_path + fname));
            }
        }
    }
    else if(v.is_array()) // handle an array of tables
    {
        for(auto& elem : v.as_array())
        {
            expand_include(elem);
        }
    }
    return;
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
