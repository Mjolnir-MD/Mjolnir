#ifndef MJOLNIR_READ_OBSERVER
#define MJOLNIR_READ_OBSERVER
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
Observer<traitsT>
read_observer(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_observer(), 0);

    const auto& files  = toml::find<toml::value>(root, "files");
    const auto& output = toml::find<toml::value>(files, "output");

    std::string output_path("./");
    if(toml::get<toml::table>(output).count("path") == 1)
    {
        output_path = toml::find<std::string>(output, "path");
        if(output_path.back() != '/') {output_path += '/';}
    }
    const auto output_prefix = toml::find<std::string>(output, "prefix");
    MJOLNIR_LOG_NOTICE("output files are `", output_path, output_prefix, ".*`");

    return Observer<traitsT>(output_path + output_prefix);
}

} // mjolnir
#endif// MJOLNIR_READ_OBSERVER
