#ifndef MJOLNIR_READ_OBSERVER
#define MJOLNIR_READ_OBSERVER
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/Observer.hpp>

namespace mjolnir
{

template<typename traitsT>
Observer<traitsT>
read_observer(const toml::Table& data)
{
    const auto& general = data.at("general").cast<toml::value_t::Table>();
    const std::string path  = toml::get<std::string>(general.at("output_path"));
    const std::string fname = toml::get<std::string>(general.at("file_name"));

    const auto& simulator = data.at("simulator").cast<toml::value_t::Table>();
    const std::size_t interval =
        toml::get<std::size_t>(simulator.at("observe_interval"));

    return Observer<traitsT>(path + fname + std::string(".xyz"),
                             path + fname + std::string(".ene"),
                             interval);
}



}
#endif// MJOLNIR_READ_OBSERVER
