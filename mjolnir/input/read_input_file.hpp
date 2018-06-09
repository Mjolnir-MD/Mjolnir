#ifndef MJOLNIR_READ_INPUT_FILE
#define MJOLNIR_READ_INPUT_FILE
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/input/read_simulator.hpp>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_parameter(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_parameter(const toml::Table& data), 0);

    typedef typename traitsT::real_type real_type;
    const auto& parameter = toml_value_at(data, "parameters", "<root>"
            ).template cast<toml::value_t::Table>();

    physics<real_type>::kB = toml::get<real_type>(
            toml_value_at(parameter, "kB", "[parameters]"));
    physics<real_type>::NA = toml::get<real_type>(
            toml_value_at(parameter, "NA", "[parameters]"));
    physics<real_type>::e  = toml::get<real_type>(
            toml_value_at(parameter, "e",  "[parameters]"));
    physics<real_type>::vacuum_permittivity = toml::get<real_type>(
            toml_value_at(parameter, "ε0", "[parameters]"));

    MJOLNIR_LOG_INFO("kB = ", physics<real_type>::kB);
    MJOLNIR_LOG_INFO("NA = ", physics<real_type>::NA);
    MJOLNIR_LOG_INFO("e  = ", physics<real_type>::e);
    MJOLNIR_LOG_INFO("ε0 = ", physics<real_type>::vacuum_permittivity);

    return read_simulator<traitsT>(data);
}

template<typename realT>
std::unique_ptr<SimulatorBase>
read_boundary(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_boundary(const toml::Table& data), 0);

    const auto& general = toml_value_at(data, "general", "<root>"
            ).template cast<toml::value_t::Table>();
    const auto boundary = toml::get<std::string>(
            toml_value_at(general, "boundary", "[general]"));

    if(boundary == "Unlimited")
    {
        MJOLNIR_LOG_INFO("boundary is UnlimitedBoundary");
        return read_parameter<
            SimulatorTraitsBase<realT, UnlimitedBoundary>>(data);
    }
    else if(boundary == "PeriodicCube")
    {
        MJOLNIR_LOG_INFO("boundary is CubicPeriodicBoudanry");
        return read_parameter<
            SimulatorTraitsBase<realT, CubicPeriodicBoundary>>(data);
    }
    else
    {
        throw std::runtime_error(
                "invalid boundary setting (Unlimited|PeriodicCube): " + boundary);
    }
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::Table& data)
{
    const auto& general = toml_value_at(data, "general", "<root>"
            ).template cast<toml::value_t::Table>();
    const auto prec = toml::get<std::string>(
            toml_value_at(general, "precision", "[general]"));

    MJOLNIR_SET_DEFAULT_LOGGER(std::string(toml::get<std::string>(
            toml_value_at(general, "file_name", "[general]")) + ".log"));
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_precision(const toml::Table& data), 0);

    if(prec == "double")
    {
        MJOLNIR_LOG_INFO("precision is double");
        return read_boundary<double>(data);
    }
    else if(prec == "float")
    {
        MJOLNIR_LOG_INFO("precision is float");
        return read_boundary<float>(data);
    }
    else
    {
        throw std::runtime_error(
                "invalid precision setting (double|float): " + prec);
    }
}

inline std::unique_ptr<SimulatorBase>
read_input_file(const std::string& filename)
{
    const auto data = toml::parse(filename);
    return read_precision(data); // read all the settings recursively...
}


}// mjolnir
#endif// MJOLNIR_READ_INPUT_FILE
