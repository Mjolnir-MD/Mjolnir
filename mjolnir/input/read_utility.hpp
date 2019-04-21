#ifndef MJOLNIR_INPUT_READ_UTILITY_HPP
#define MJOLNIR_INPUT_READ_UTILITY_HPP
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <extlib/toml/toml.hpp>
#include <type_traits>
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace mjolnir
{

// this function reads array of structs from a toml value like the following.
// ```toml
// v = [
//    {index = 0, a = 1.0, b = 100},
//    {index = 1, a = 1.0, b = 100},
//    ...
//    {index =99, a = 1.0, b = 100},
// ]
// ```
// or 
// ```toml
// v.size    = 100                    # length
// v.default = {a = 1.0, b = 100}     # default value
// v.values  = [
//     {index = 10, a = 1.0, b = 100} # overwrite special values
// ]
// ```
//
// In order to read a table like `{a = 1.0, b = 100}` as a struct like
// `std::pair<double, double>`, we need to pass a function that converts
// `toml::value` to `std::pair<double, double>`. So the function receives
// `reader(const toml::value&) -> T`.
//
// Thus the reader MUST ignore a key named "index".
//
// Since we need to convert several kinds of tables that has different keys
// into the same struct like `std::pair<double, double>`, `toml::from` is
// not helpful because it can only be implemented once for one struct.
template<typename T, typename F>
std::vector<T> read_array(const toml::value& v, F&& reader)
{
    static_assert(std::is_default_constructible<T>::value, "");
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // This function is for logging. It always formats toml value as inline.
    const auto format_inline = [](const toml::value& v) -> std::string {
        return toml::format(v,
            /* width      = */ std::numeric_limits<std::size_t>::max(),
            /* float prec = */ 6, /* force_inline = */true);
    };

    if(v.is_array())
    {
        const auto& arr = toml::get<toml::array>(v);
        if(arr.empty())
        {
            return std::vector<T>{};
        }

        std::vector<T> vec(arr.size());
        MJOLNIR_LOG_INFO(arr.size(), " parameters found.");

        if(arr.front().is_table())
        {
            std::string msg;
            bool maybe_exhaustive = true;
            std::vector<std::size_t> idxs; idxs.reserve(arr.size());

            for(const toml::value& item : arr)
            {
                const auto i = toml::find<std::size_t>(item, "index");
                if(i >= vec.size())
                {
                    // the index exceeds the size of the array. it cannot be
                    // exhaustive. break the loop and throw an error.
                    maybe_exhaustive = false;
                    msg = toml::format_error("the index exceeds the number of "
                                             "parameters", item, "here");
                    break;
                }
                vec.at(i) = reader(item);
                idxs.push_back(i);

                MJOLNIR_LOG_INFO(i, "-th parameter, ", format_inline(item),
                                 ", has been read.");
            }

            if(maybe_exhaustive)
            {
                // check all the elements are correctly initialized
                if(!std::is_sorted(idxs.begin(), idxs.end()))
                {
                    std::sort(idxs.begin(), idxs.end());
                }
                const auto uniqued = std::unique(idxs.begin(), idxs.end());
                if(uniqued == idxs.end())
                {
                    // all the values are initialized. the values are okay.
                    return vec;
                }
                // otherwise, not okay. raise an error.
                std::ostringstream oss;
                oss << "indexes ";
                for(auto iter = uniqued; iter!=idxs.end(); ++iter)
                {
                    oss << *iter << ' ';
                }
                oss << "are overlapping.";
                msg = oss.str();
            }
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_array: values are not exhaustive.",
                v, "index may overlap or exceeds an range.", {msg}));
        }
        // otherwise, format error. throw an exception.
    }
    else if(v.is_table())
    {
        // v.size    = 100                    # length
        // v.default = {a = 1.0, b = 100}     # default value
        // v.values  = [
        //     {index = 10, a = 1.0, b = 100} # overwrite special values
        // ]
        const auto& tab = toml::get<toml::table>(v);
        if(tab.count("size") == 1 && tab.count("default") == 1)
        {
            const auto len  = toml::find<std::size_t>(v, "size");
            const auto dflt = reader(toml::find<toml::value>(v, "default"));
            std::vector<T> vec(len, dflt);

            MJOLNIR_LOG_INFO("number of parameters are ", len);
            MJOLNIR_LOG_INFO("values of default parameter is ",
                             format_inline(toml::find<toml::value>(v, "default")));

            if(tab.count("values") == 1)
            {
                for(const auto& sp_val : toml::find<toml::array>(v, "values"))
                {
                    const auto idx = toml::find<std::size_t>(sp_val, "index");
                    if(idx >= len)
                    {
                        throw_exception<std::runtime_error>(toml::format_error(
                            "[error] mjolnir::read_array: the index exceeds the"
                            " length of array", sp_val, "here"));
                    }
                    vec.at(idx) = reader(sp_val);

                    MJOLNIR_LOG_INFO("overwrite ", idx, "-th parameter by ",
                                     format_inline(sp_val));
                }
            }
            return vec;
        }
        // otherwise (if no `size` and `default` exist), raise format error.
    }
    throw_exception<std::runtime_error>(toml::format_error("[error] "
        "mjolnir::read_array: invalid format for an array of parameters",
        v, "here", {"expected array of parameters. the format is:",
        "key = [ # write all parameters explicitly",
        "    {index=0, a=1.0, b=100},",
        "    ...",
        "]",
        "or",
        "key.size    = 100",
        "key.default = {a=1.0, b=100}",
        "key.values  = [",
        "    {index = 10, a = 2.0, b=10}, # overwrite parameter for i = 10",
        "    ...",
        "]"
        }));
}


} // mjolnir
#endif// MJOLNIR_INPUT_READ_UTILITY_HPP
