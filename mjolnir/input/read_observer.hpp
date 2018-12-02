#ifndef MJOLNIR_READ_OBSERVER
#define MJOLNIR_READ_OBSERVER
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
Observer<traitsT>
read_observer(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_observer(), 0);

    const auto& files =
        get_toml_value<toml::Table>(data, "files", "<root>");
    const auto& output =
        get_toml_value<toml::Table>(files, "output", "[files]");

    std::string output_path("./");
    if(output.count("path") == 1)
    {
        output_path = get_toml_value<std::string>(output, "path", "[files.output]");
    }
    if(output_path.back() != '/') {output_path += '/';}

    const auto output_prefix = get_toml_value<std::string>(
            output, "prefix", "[files.output]");

    MJOLNIR_LOG_INFO("path   = ", output_path);
    MJOLNIR_LOG_INFO("prefix = ", output_prefix);

    MJOLNIR_LOG_NOTICE("output files are `", output_path, output_prefix, ".*`");

    return Observer<traitsT>(output_path + output_prefix);
}

} // mjolnir
#endif// MJOLNIR_READ_OBSERVER
