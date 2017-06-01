#ifndef MJOLNIR_READ_INPUT_FILE
#define MJOLNIR_READ_INPUT_FILE
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include "read_simulator.hpp"
#include <memory>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_parameter(const toml::Table& data)
{
    typedef typename traitsT::real_type real_type;
    const auto& parameter =
        data.at("parameter").template cast<toml::value_t::Table>();

    physics<real_type>::kB = toml::get<real_type>(parameter.at("kB"));
    physics<real_type>::NA = toml::get<real_type>(parameter.at("NA"));
    physics<real_type>::e  = toml::get<real_type>(parameter.at("e"));
    physics<real_type>::vacuum_permittivity =
        toml::get<real_type>(parameter.at("vacuum_permittivity"));

    return read_simulator<traitsT>(data);
}

template<typename realT>
std::unique_ptr<SimulatorBase>
read_boundary(const toml::Table& data)
{
    const auto& general = data.at("general").template cast<toml::value_t::Table>();
    std::string boundary = toml::get<std::string>(general.at("boundary"));
    if(boundary == "Unlimited")
    {
        return read_parameter<
            SimulatorTraitsBase<realT, UnlimitedBoundary>>(data);
    }
    else if(boundary == "PeriodicCube")
    {
        return read_parameter<
            SimulatorTraitsBase<realT, CubicPeriodicBoundary>>(data);
    }
    else
    {
        throw std::runtime_error(
                "invalid boundary setting (double|float): " + boundary);
    }
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::Table& data)
{
    const auto& general = data.at("general").template cast<toml::value_t::Table>();
    std::string prec = toml::get<std::string>(general.at("precision"));
    if(prec == "double")
    {
        return read_boundary<double>(data);
    }
    else if(prec == "float")
    {
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
