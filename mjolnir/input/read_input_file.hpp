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
        MJOLNIR_LOG_NOTICE("execute on single core");
        return read_units<SimulatorTraits<realT, boundaryT>>(root, simulator);
    }

    const auto parallelism = toml::find(simulator, "parallelism");
    if(parallelism.is_string() && parallelism.as_string() == "sequencial")
    {
        MJOLNIR_LOG_NOTICE("execute on single core");
        return read_units<SimulatorTraits<realT, boundaryT>>(root, simulator);
    }
    else if(parallelism.is_string() &&
           (parallelism.as_string() == "openmp" ||
            parallelism.as_string() == "OpenMP"))
    {
#ifdef MJOLNIR_WITH_OPENMP
        MJOLNIR_LOG_NOTICE("execute on ", omp_get_max_threads() ," cores with openmp");
        return read_units<OpenMPSimulatorTraits<realT, boundaryT>>(root, simulator);
#else
        MJOLNIR_LOG_WARN("OpenMP flag is set, but OpenMP is not enabled when building.");
        MJOLNIR_LOG_WARN("Cannot use OpenMP, running with single core.");
        return read_units<SimulatorTraits<realT, boundaryT>>(root, simulator);
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

#ifndef MJOLNIR_WITHOUT_UNLIMITED_BOUNDARY
    if(boundary == "Unlimited")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is Unlimited");
        return read_parallelism<realT, UnlimitedBoundary>(root, simulator);
    }
#endif

#ifndef MJOLNIR_WITHOUT_PERIODIC_BOUNDARY
    if(boundary == "Periodic")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is Periodic. "
                           "The shape is cuboid.");
        return read_parallelism<realT, CuboidalPeriodicBoundary>(root, simulator);
    }
    else if(boundary == "PeriodicCuboid")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is PeriodicCuboid");
        return read_parallelism<realT, CuboidalPeriodicBoundary>(root, simulator);
    }
#endif

    throw_exception<std::runtime_error>(toml::format_error("[error] "
        "mjolnir::read_boundary: invalid boundary",
        toml::find(simulator, "boundary_type"), "here", {
        "expected value is one of the following."
#ifndef MJOLNIR_WITHOUT_UNLIMITED_BOUNDARY
        , "- \"Unlimited\"     : no boundary condition. infinite space"
#endif
#ifndef MJOLNIR_WITHOUT_PERIODIC_BOUNDARY
        , "- \"Periodic\"      : periodic boundary. Assuming cuboidal shape."
        , "- \"PeriodicCuboid\": periodic boundary with cuboidal shape"
#endif
        }));
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto prec = toml::find<std::string>(simulator, "precision");

#ifndef MJOLNIR_WITHOUT_DOUBLE_PRECISION
    if(prec == "double")
    {
        MJOLNIR_LOG_NOTICE("precision is double");
        return read_boundary<double>(root, simulator);
    }
#endif

#ifndef MJOLNIR_WITHOUT_SINGLE_PRECISION
    if(prec == "float")
    {
        MJOLNIR_LOG_NOTICE("precision is float");
        return read_boundary<float>(root, simulator);
    }
#endif

    throw_exception<std::runtime_error>(toml::format_error("[error] "
        "mjolnir::read_precision: invalid precision",
        toml::find(simulator, "precision"), "here", {
        "expected value is one of the following."
#ifndef MJOLNIR_WITHOUT_DOUBLE_PRECISION
        , "- \"double\": 64 bit floating-point"
#endif
#ifndef MJOLNIR_WITHOUT_SINGLE_PRECISION
        , "- \"float\" : 32 bit floating-point"
#endif
        }));
}

inline std::unique_ptr<SimulatorBase>
read_input_file(const std::string& filename)
{
    // here, logger name is not given yet. output status directory on console.
    std::cerr << "-- reading and parsing toml file `" << filename << "` ... ";
    auto root = toml::parse(filename);
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

    MJOLNIR_LOG_NOTICE("mjolnir version ", MJOLNIR_VERSION);
    MJOLNIR_LOG_NOTICE("compiled using ",  MJOLNIR_COMPILER_VERSION);

    // Check top-level toml-values. Since it uses logger to warn,
    // we need to call it after `MJOLNIR_SET_DEFAULT_LOGGER(logger_name)`.
    check_keys_available(root, {"files"_s, "units"_s, "simulator"_s,
                                "systems"_s, "forcefields"_s});

    // load the input path in the [files] table and set it globally
    read_input_path(root);

    MJOLNIR_LOG_NOTICE("expanding include files ...");
    expand_include(root);
    MJOLNIR_LOG_NOTICE("done.");

    // the most of important flags are defined in [simulator], like
    // `precision = "float"`, `boundary_type = "Unlimited"`.
    // Those values should be read before others.
    const auto simulator = read_table_from_file(
            toml::find(root, "simulator"), "simulator");

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
