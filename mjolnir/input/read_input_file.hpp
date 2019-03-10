#ifndef MJOLNIR_READ_INPUT_FILE
#define MJOLNIR_READ_INPUT_FILE
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_units.hpp>
#include <mjolnir/input/read_files_table.hpp>
#include <memory>

namespace mjolnir
{

template<typename realT>
std::unique_ptr<SimulatorBase>
read_boundary(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // [simulator] can be provided in a different file. in that case, the table
    // has `file_name` field. In that case, file.input.path is also needed to
    // determine the location of the file.
    toml::value boundary;
    const auto& simulator = toml::find(root, "simulator");
    if(toml::get<toml::table>(simulator).count("file_name") == 1)
    {
        const auto filename   = toml::find<std::string>(simulator, "file_name");
        const auto input_path = read_input_path(root);
        const auto simfile    = toml::parse(input_path + filename);
        if(simfile.count("simulator") != 1)
        {
            throw_exception<std::out_of_range>("[error] mjolnir::read_simulator: "
                "table [simulator] not found in the toml file\n --> ",
                input_path, filename, "\n | the file should define [simulator] "
                "table and define values in it.");
        }
        boundary = toml::find(simfile.at("simulator"), "boundary_type");
    }
    else
    {
        boundary = toml::find(simulator, "boundary_type");
    }

    if(toml::get<std::string>(boundary) == "Unlimited")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is Unlimited");
        return read_units<SimulatorTraits<realT, UnlimitedBoundary>>(root);
    }
    else if(toml::get<std::string>(boundary) == "PeriodicCuboid")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is CuboidalPeriodic");
        return read_units<SimulatorTraits<realT, CuboidalPeriodicBoundary>>(root);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_boundary: invalid boundary", boundary, "here", {
            "- \"Unlimited\"     : no boundary condition. infinite space",
            "- \"PeriodicCuboid\": periodic boundary with cuboidal shape"
            }));
    }
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    // [simulator] can be provided in a different file. in that case, the table
    // has `file_name` field. In that case, input_path is also needed to
    // determine the location of the file.
    toml::value prec;

    const auto& simulator = toml::find(root, "simulator");
    if(toml::get<toml::table>(simulator).count("file_name") == 1)
    {
        const auto filename   = toml::find<std::string>(simulator, "file_name");
        const auto input_path = read_input_path(root);
        const auto simfile    = toml::parse(input_path + filename);
        if(simfile.count("simulator") != 1)
        {
            throw_exception<std::out_of_range>("[error] mjolnir::read_simulator: "
                "table [simulator] not found in the toml file\n --> ",
                input_path, filename, "\n | the file should define [simulator] "
                "table and define values in it.");
        }
        prec = toml::find(simfile.at("simulator"), "precision");
    }
    else
    {
        prec = toml::find(simulator, "precision");
    }

    if(toml::get<std::string>(prec) == "double")
    {
        MJOLNIR_LOG_NOTICE("precision is double");
        return read_boundary<double>(root);
    }
    else if(toml::get<std::string>(prec) == "float")
    {
        MJOLNIR_LOG_NOTICE("precision is float");
        return read_boundary<float>(root);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_precision: invalid precision", prec, "here", {
            "expected value is one of the following.",
            "- \"double\": 64 bit floating-point",
            "- \"float\" : 32 bit floating-point"
            }));
    }
}

inline std::unique_ptr<SimulatorBase>
read_input_file(const std::string& filename)
{
    // here, logger name is not given yet. output status directory on console.
    std::cerr << "-- reading and parsing toml file `" << filename << "` ... ";
    const auto root = toml::parse(filename);
    std::cerr << " successfully parsed." << std::endl;

    // initializing logger by using output_path and output_prefix ...
    const auto& files    = toml::find<toml::value>(root,   "files");
    const auto& output   = toml::find<toml::value>(files,  "output");
    const auto  out_path = toml::find<std::string>(output, "path");

    // XXX:  Here, this code assumes POSIX. it does not support windows.
    // TODO: Consider using Boost.filesystem to manage path and files
    //       in more elegant and powerful way? After switching C++17,
    //       we can re-write that with <filesystem>.
    const auto logger_name = out_path +
        toml::find<std::string>(output, "prefix") + ".log";

    MJOLNIR_SET_DEFAULT_LOGGER(logger_name);
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_NOTICE("the log file is `", logger_name, '`');

    return read_precision(root); // read all the settings recursively...
}

}// mjolnir
#endif// MJOLNIR_READ_INPUT_FILE
