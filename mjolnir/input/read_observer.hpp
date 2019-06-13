#ifndef MJOLNIR_INPUT_READ_OBSERVER_HPP
#define MJOLNIR_INPUT_READ_OBSERVER_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/core/ObserverContainer.hpp>
#include <mjolnir/core/EnergyObserver.hpp>
#include <mjolnir/core/XYZObserver.hpp>
#include <mjolnir/core/DCDObserver.hpp>
#include <mjolnir/core/TRRObserver.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/make_unique.hpp>

namespace mjolnir
{

template<typename traitsT>
void add_observer(ObserverContainer<traitsT>& observers,
                  const toml::value& format, const std::string& file_prefix)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // To show the informative error message, here it uses toml::value that
    // contains file location. But this function assumes that the `format`
    // contains `toml::string`.
    assert(format.is(toml::value_t::String));

    if(format == "xyz")
    {
        using observer_type = XYZObserver<traitsT>;
        MJOLNIR_LOG_NOTICE("output xyz format.");
        observers.push_back(make_unique<observer_type>(file_prefix));
        return;
    }
    else if(format == "dcd")
    {
        using observer_type = DCDObserver<traitsT>;
        MJOLNIR_LOG_NOTICE("output dcd format.");
        observers.push_back(make_unique<observer_type>(file_prefix));
        return;
    }
    else if(format == "trr")
    {
        using observer_type = TRRObserver<traitsT>;
        MJOLNIR_LOG_NOTICE("output trr format.");
        observers.push_back(make_unique<observer_type>(file_prefix));
        return;
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_observer: output format not supported", format,
            "here", {
                "expected one of the following.",
                "- \"xyz\": the simplest ascii format.",
                "- \"dcd\": binary format which contains normally positions."
                "- \"trr\": binary format which contains positions, velocities, and forces."
            }));
    }
    return ;
}


template<typename traitsT>
ObserverContainer<traitsT>
read_observer(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto& files  = toml::find<toml::value>(root, "files");
    const auto& output = toml::find<toml::value>(files, "output");

    const auto progress_bar_enabled = toml::expect<bool>(output, "progress_bar")
                                                        .unwrap_or(true);

    const auto output_path   = read_output_path(root);
    const auto output_prefix = toml::find<std::string>(output, "prefix");
    MJOLNIR_LOG_NOTICE("output file prefix is `", output_path, output_prefix, '`');

    const std::string file_prefix = output_path + output_prefix;

    // ------------------------------------------------------------------------

    ObserverContainer<traitsT> observers(progress_bar_enabled);

    // Energy is always written to "prefix.ene".
    observers.push_back(make_unique<EnergyObserver<traitsT>>(file_prefix));

    const auto& format = toml::find(output, "format");

    if(format.is(toml::value_t::String))
    {
        add_observer(observers, format, file_prefix);
    }
    else if(format.is(toml::value_t::Array))
    {
        const auto fmts = toml::get<toml::array>(format);
        for(const auto& fmt : fmts)
        {
            add_observer(observers, fmt, file_prefix);
        }
    }
    return observers;
}

} // mjolnir
#endif// MJOLNIR_READ_OBSERVER_HPP
