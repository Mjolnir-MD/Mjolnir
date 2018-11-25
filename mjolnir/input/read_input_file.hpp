#ifndef MJOLNIR_READ_INPUT_FILE
#define MJOLNIR_READ_INPUT_FILE
#include <extlib/toml/toml/toml.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_units.hpp>
#include <memory>

namespace mjolnir
{

template<typename realT>
std::unique_ptr<SimulatorBase>
read_boundary(const toml::table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_boundary(const toml::table& data), 0);

    // [simulator] can be provided in a different file. in that case, the table
    // has `file_name` field. In that case, input_path is also needed to
    // determine the location of the file.

    std::string boundary;

    const auto& simulator =
        get_toml_value<toml::table>(data, "simulator", "<root>");

    if(simulator.count("file_name") == 1)
    {
        const auto filename = toml::get<std::string>(simulator.at("file_name"));
        std::string input_path; //default: empty

        const auto& files = get_toml_value<toml::table>(data, "files", "<root>");
        if(files.count("input_path") == 1)
        {
            input_path = toml::get<std::string>(files.at("input_path"));
        }

        const auto& sim = toml::parse(input_path + filename);
        boundary = get_toml_value<std::string>(
                sim, "boundary_type", "[simulator]");
    }
    else
    {
        boundary = get_toml_value<std::string>(
                simulator, "boundary_type", "[simulator]");
    }

    if(boundary == "Unlimited")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is Unlimited");
        return read_units<SimulatorTraits<realT, UnlimitedBoundary>>(data);
    }
    else if(boundary == "PeriodicCuboid")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is CuboidalPeriodic");
        return read_units<SimulatorTraits<realT, CuboidalPeriodicBoundary>>(data);
    }
    else
    {
        throw_exception<std::runtime_error>("mjolnir::read_boundary: "
            "invalid boundary: neither 'Unlimited' or 'PeriodicCuboid' -> ",
            boundary);
    }
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_precision(const toml::table& data), 0);

    // [simulator] can be provided in a different file. in that case, the table
    // has `file_name` field. In that case, input_path is also needed to
    // determine the location of the file.

    std::string prec;

    const auto& simulator =
        get_toml_value<toml::table>(data, "simulator", "<root>");

    if(simulator.count("file_name") == 1)
    {
        const auto filename = toml::get<std::string>(simulator.at("file_name"));
        std::string input_path; //default: empty

        const auto& files = get_toml_value<toml::table>(data, "files", "<root>");
        if(files.count("input_path") == 1)
        {
            input_path = toml::get<std::string>(files.at("input_path"));
        }

        const auto& sim = toml::parse(input_path + filename);
        prec = get_toml_value<std::string>(sim, "precision", "[simulator]");
    }
    else
    {
        prec = get_toml_value<std::string>(simulator, "precision", "[simulator]");
    }

    if(prec == "double")
    {
        MJOLNIR_LOG_NOTICE("precision is double");
        return read_boundary<double>(data);
    }
    else if(prec == "float")
    {
        MJOLNIR_LOG_NOTICE("precision is float");
        return read_boundary<float>(data);
    }
    else
    {
        throw_exception<std::runtime_error>("mjolnir::read_precision: "
            "invalid precision: neither 'double' or 'float' -> ", prec);
    }
}

inline std::unique_ptr<SimulatorBase>
read_input_file(const std::string& filename)
{
    // here, logger name is not given yet. output status directory on console.
    std::cerr << "-- reading and parsing toml file `" << filename << "` ... ";
    const auto data = toml::parse(filename);
    std::cerr << " successfully parsed." << std::endl;

    // initializing logger by using output_path and output_prefix ...
    const auto& files = get_toml_value<toml::table>(data, "files", "<root>");
    const auto  path  = get_toml_value<std::string>(files, "output_path", "[files]");

    // XXX:  Here, this code assumes POSIX. it does not support windows.
    // TODO: Consider using Boost.filesystem to manage path and files
    //       in more elegant and powerful way? After switching C++17,
    //       we can re-write that with <filesystem>.

    const std::string logger_name = path + get_toml_value<std::string>(
            files, "output_prefix", "[files]") + ".log";
    MJOLNIR_SET_DEFAULT_LOGGER(logger_name);
    MJOLNIR_GET_DEFAULT_LOGGER();

    MJOLNIR_LOG_NOTICE("the log file is `", logger_name, '`');
    MJOLNIR_SCOPE(read_input_file(const toml::table& data), 0);

    return read_precision(data); // read all the settings recursively...
}

}// mjolnir
#endif// MJOLNIR_READ_INPUT_FILE
