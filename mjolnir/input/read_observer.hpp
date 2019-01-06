#ifndef MJOLNIR_READ_OBSERVER_HPP
#define MJOLNIR_READ_OBSERVER_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ObserverBase.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/make_unique.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<ObserverBase<traitsT>>
read_observer(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_observer(), 0);

    const auto& files  = toml::find<toml::value>(root, "files");
    const auto& output = toml::find<toml::value>(files, "output");

    const auto progress_bar = toml::expect<bool>(output, "progress_bar")
                                                .unwrap_or(true);

    std::string output_path("./");
    if(toml::get<toml::table>(output).count("path") == 1)
    {
        output_path = toml::find<std::string>(output, "path");
        if(output_path.back() != '/') {output_path += '/';}
    }
    const auto output_prefix = toml::find<std::string>(output, "prefix");
    MJOLNIR_LOG_NOTICE("output files are `", output_path, output_prefix, ".*`");

    const std::string file_prefix = output_path + output_prefix;

    const auto& format = toml::find<std::string>(output, "format");
    if(format == "xyz")
    {
        using observer_type = Observer<traitsT>;
        MJOLNIR_LOG_NOTICE("output format is xyz.");
        return make_unique<observer_type>(file_prefix, progress_bar);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_observer: output format not supported",
            toml::find(output, "format"), "here", {
            "expected one of the following.",
            "- \"xyz\": the simplest ascii format",
            }));
    }
}

} // mjolnir
#endif// MJOLNIR_READ_OBSERVER_HPP
