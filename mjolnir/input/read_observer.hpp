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

    std::string path_ =
        get_toml_value<std::string>(files, "output_path", "[files]");
    if(path_.back() != '/') {path_ += '/';} //XXX assuming posix
    const auto path(path_);

    const std::string prefix = get_toml_value<std::string>(
            files, "output_prefix", "[files]");
    MJOLNIR_LOG_INFO("path   = ", path);
    MJOLNIR_LOG_INFO("prefix = ", prefix);

    MJOLNIR_LOG_NOTICE("output files are `", path, prefix, ".*`");

    return Observer<traitsT>(path + prefix);
}



}
#endif// MJOLNIR_READ_OBSERVER
