#ifndef MJOLNIR_READ_OBSERVER
#define MJOLNIR_READ_OBSERVER
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/util/get_toml_value.hpp>

namespace mjolnir
{

template<typename traitsT>
Observer<traitsT>
read_observer(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_observer(), 0);

    const auto& general = toml_value_at(data, "general", "<root>"
            ).cast<toml::value_t::Table>();
    std::string path = toml::get<std::string>(
            toml_value_at(general, "output_path", "[general]"));
    //XXX assuming posix
    if(path.back() != '/') {path += '/';}

    const std::string prefix = toml::get<std::string>(
            toml_value_at(general, "output_prefix", "[general]"));
    MJOLNIR_LOG_INFO("path   = ", path);
    MJOLNIR_LOG_INFO("prefix = ", prefix);

    const auto& simulator = toml_value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    const std::size_t interval = toml::get<std::size_t>(
            toml_value_at(simulator, "save_step", "[simulator]"));
    MJOLNIR_LOG_INFO("interval = ", interval);

    return Observer<traitsT>(path + prefix, interval);
}



}
#endif// MJOLNIR_READ_OBSERVER
