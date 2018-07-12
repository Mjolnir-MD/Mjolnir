#ifndef MJOLNIR_READ_INPUT_FILE
#define MJOLNIR_READ_INPUT_FILE
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
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

    const auto& parameter =
        get_toml_value<toml::Table>(data, "parameters", "<root>");

    physics<real_type>::kB = get_toml_value<real_type>(parameter, "kB", "[parameters]");
    physics<real_type>::NA = get_toml_value<real_type>(parameter, "NA", "[parameters]");
    physics<real_type>::e  = get_toml_value<real_type>(parameter, "e",  "[parameters]");
    physics<real_type>::vacuum_permittivity =
        get_toml_value<real_type>(parameter, "ε0", "[parameters]");

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

    const auto& general = get_toml_value<toml::Table>(data, "general", "<root>");
    const auto boundary =
        get_toml_value<std::string>(general, "boundary", "[general]");

    if(boundary == "Unlimited")
    {
        MJOLNIR_LOG_INFO("boundary is UnlimitedBoundary");
        return read_parameter<
            SimulatorTraitsBase<realT, UnlimitedBoundary>>(data);
    }
    else if(boundary == "PeriodicCuboid")
    {
        MJOLNIR_LOG_INFO("boundary is CuboidalPeriodicBoudanry");
        return read_parameter<
            SimulatorTraitsBase<realT, CuboidalPeriodicBoundary>>(data);
    }
    else
    {
        throw std::runtime_error(
            "invalid boundary setting (Unlimited|PeriodicCuboid): " + boundary);
    }
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_precision(const toml::Table& data), 0);

    const auto& general = get_toml_value<toml::Table>(data, "general", "<root>");
    const auto prec =
        get_toml_value<std::string>(general, "precision", "[general]");

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

    // setting logger ...
    const auto& general = get_toml_value<toml::Table>(data, "general", "<root>");
    std::string path = get_toml_value<std::string>(general, "output_path", "[general]");
    if(path.back() != '/') {path += '/';/*XXX assuming posix */}

    const std::string logger_name = path + get_toml_value<std::string>(
            general, "output_prefix", "[general]") + ".log";
    MJOLNIR_SET_DEFAULT_LOGGER(logger_name);

    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_input_file(const toml::Table& data), 0);

    return read_precision(data); // read all the settings recursively...
}

}// mjolnir
#endif// MJOLNIR_READ_INPUT_FILE
