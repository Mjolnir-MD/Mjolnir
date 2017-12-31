#ifndef MJOLNIR_READ_OBSERVER
#define MJOLNIR_READ_OBSERVER
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/input/get_toml_value.hpp>

namespace mjolnir
{

template<typename traitsT>
Observer<traitsT>
read_observer(const toml::Table& data)
{
    const auto& general = toml_value_at(data, "general", "<root>"
            ).cast<toml::value_t::Table>();
    const std::string path  = toml::get<std::string>(
            toml_value_at(general, "output_path", "[general]"));
    const std::string fname = toml::get<std::string>(
            toml_value_at(general, "file_name", "[general]"));

    const auto& simulator = toml_value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    const std::size_t interval = toml::get<std::size_t>(
            toml_value_at(simulator, "save_step", "[simulator]"));

    return Observer<traitsT>(path + fname + std::string(".xyz"),
                             path + fname + std::string(".ene"),
                             interval);
}



}
#endif// MJOLNIR_READ_OBSERVER
