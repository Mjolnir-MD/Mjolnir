#ifndef MJOLNIR_INPUT_READ_INPUT_FILE_HPP
#define MJOLNIR_INPUT_READ_INPUT_FILE_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/input/read_units.hpp>
#include <mjolnir/input/read_path.hpp>

#ifdef MJOLNIR_WITH_OPENMP
#include <mjolnir/omp/omp.hpp>
#endif

#include <memory>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
std::unique_ptr<SimulatorBase>
read_parallelism(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if(simulator.as_table().count("parallelism") == 0)
    {
        return read_units<SimulatorTraits<realT, boundaryT>>(root);
    }

    const auto parallelism = toml::find(simulator, "parallelism");
    if(parallelism.is_string() && parallelism.as_string() == "sequencial")
    {
        MJOLNIR_LOG_NOTICE("execute on single core");
        return read_units<SimulatorTraits<realT, boundaryT>>(root);
    }
    else if(parallelism.is_string() &&
           (parallelism.as_string() == "openmp" ||
            parallelism.as_string() == "OpenMP"))
    {
#ifdef MJOLNIR_WITH_OPENMP
        MJOLNIR_LOG_NOTICE("execute on ", omp_get_max_threads() ," cores with openmp");
        return read_units<OpenMPSimulatorTraits<realT, boundaryT>>(root);
#else
        MJOLNIR_LOG_WARN("OpenMP flag is set, but OpenMP is not enabled when building.");
        MJOLNIR_LOG_WARN("Cannot use OpenMP, running with single core.");
        return read_units<SimulatorTraits<realT, boundaryT>>(root);
#endif
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_parallelism: invalid variable ",
            toml::find(simulator, "parallelism"), "here", {
            "- \"sequencial\": run with only 1 core (default)",
            "- \"openmp\"    : use openmp to parallelize."
            }));
    }
}

template<typename realT>
std::unique_ptr<SimulatorBase>
read_boundary(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto boundary = toml::find<std::string>(simulator, "boundary_type");
    if(boundary == "Unlimited")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is Unlimited");
        return read_parallelism<realT, UnlimitedBoundary>(root, simulator);
    }
    else if(boundary == "PeriodicCuboid")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is CuboidalPeriodic");
        return read_parallelism<realT, CuboidalPeriodicBoundary>(root, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_boundary: invalid boundary",
            toml::find(simulator, "boundary_type"), "here", {
            "- \"Unlimited\"     : no boundary condition. infinite space",
            "- \"PeriodicCuboid\": periodic boundary with cuboidal shape"
            }));
    }
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto prec = toml::find<std::string>(simulator, "precision");
    if(prec == "double")
    {
        MJOLNIR_LOG_NOTICE("precision is double");
        return read_boundary<double>(root, simulator);
    }
    else if(prec == "float")
    {
        MJOLNIR_LOG_NOTICE("precision is float");
        return read_boundary<float>(root, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_precision: invalid precision",
            toml::find(simulator, "precition"), "here", {
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
    const auto& output   = toml::find(root, "files", "output");
    const auto  out_path = read_output_path(root);

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

    // Check top-level toml-values. Since it uses logger to warn,
    // we need to call it after `MJOLNIR_SET_DEFAULT_LOGGER(logger_name)`.
    check_keys_available(root, {"files"_s, "units"_s, "simulator"_s,
                                "systems"_s, "forcefields"_s});

    // the most of important flags are defined in [simulator], like
    // `precision = "float"`, `boundary_type = "Unlimited"`.
    // Those values should be read before others.
    // Thus first read [simulator] here and pass it to the latter functions.

    const auto& simulator = toml::find(root, "simulator");
    if(simulator.as_table().count("file_name") == 1)
    {
        // [simulator] can be provided in a different file. In that case, the
        // table has `file_name` field. In that case, input_path is also needed
        // to determine the location of the file.
        const auto filename   = toml::find<std::string>(simulator, "file_name");
        const auto input_path = read_input_path(root);
        const auto simfile    = toml::parse(input_path + filename);
        if(simfile.as_table().count("simulator") != 1)
        {
            throw_exception<std::out_of_range>("[error] mjolnir::read_simulator: "
                "table [simulator] not found in the toml file\n --> ",
                input_path, filename, "\n | the file should define [simulator] "
                "table and define values in it.");
        }
        return read_precision(root, simfile.as_table().at("simulator"));
    }
    return read_precision(root, simulator);
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<SimulatorBase> read_parallelism<double, UnlimitedBoundary       >(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_parallelism<float , UnlimitedBoundary       >(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_parallelism<double, CuboidalPeriodicBoundary>(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_parallelism<float , CuboidalPeriodicBoundary>(const toml::value& root, const toml::value& simulator);

extern template std::unique_ptr<SimulatorBase> read_boundary<double>(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_boundary<float >(const toml::value& root, const toml::value& simulator);
#endif

}// mjolnir
#endif// MJOLNIR_READ_INPUT_FILE
