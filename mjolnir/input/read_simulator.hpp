#ifndef MJOLNIR_INPUT_READ_SIMULATOR_HPP
#define MJOLNIR_INPUT_READ_SIMULATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/MolecularDynamicsSimulator.hpp>
#include <mjolnir/core/SteepestDescentSimulator.hpp>
#include <mjolnir/core/SimulatedAnnealingSimulator.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_system.hpp>
#include <mjolnir/input/read_forcefield.hpp>
#include <mjolnir/input/read_integrator.hpp>
#include <mjolnir/input/read_observer.hpp>
#include <mjolnir/input/read_path.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_molecular_dynamics_simulator(
        const toml::table& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto& integrator     = toml::find(simulator, "integrator");
    const auto integrator_type = toml::find<std::string>(integrator, "type");
    const auto tstep      = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep      = toml::find<std::size_t>(simulator, "save_step");
    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    if(integrator_type == "VelocityVerlet")
    {
        MJOLNIR_LOG_NOTICE("Integrator is VelocityVerlet.");
        using integrator_t = VelocityVerletIntegrator<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;

        return make_unique<simulator_t>(
                tstep, sstep,
                read_system<traitsT>(root, 0),
                read_forcefield<traitsT>(root, 0),
                read_velocity_verlet_integrator<traitsT>(simulator),
                read_observer<traitsT>(root));
    }
    else if(integrator_type == "UnderdampedLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is Underdamped Langevin.");
        using integrator_t = UnderdampedLangevinIntegrator<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;

        return make_unique<simulator_t>(
                tstep, sstep,
                read_system<traitsT>(root, 0),
                read_forcefield<traitsT>(root, 0),
                read_underdamped_langevin_integrator<traitsT>(simulator),
                read_observer<traitsT>(root));
    }
    else if(integrator_type == "BAOABLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is BAOAB Langevin.");
        using integrator_t = BAOABLangevinIntegrator<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;

        return make_unique<simulator_t>(
                tstep, sstep,
                read_system<traitsT>(root, 0),
                read_forcefield<traitsT>(root, 0),
                read_BAOAB_langevin_integrator<traitsT>(simulator),
                read_observer<traitsT>(root));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_molecular_dynamics_simulator: invalid integrator: ",
            toml::find(integrator, "type"), "here", {
            "expected value is one of the following.",
            "- \"VelocityVerlet\"     : simple and standard Velocity Verlet integrator.",
            "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                      " based on the Velocity Verlet"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_steepest_descent_simulator(
        const toml::table& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using simulator_type = SteepestDescentSimulator<traitsT>;

    const auto step_lim  = toml::find<std::size_t>(simulator, "step_limit");
    const auto save_step = toml::find<std::size_t>(simulator, "save_step");
    const auto delta     = toml::find<real_type  >(simulator, "delta");
    const auto threshold = toml::find<real_type  >(simulator, "threshold");

    MJOLNIR_LOG_NOTICE("step_limit is ", step_lim);
    MJOLNIR_LOG_NOTICE("save_step  is ", save_step);
    MJOLNIR_LOG_NOTICE("delta      is ", delta);
    MJOLNIR_LOG_NOTICE("threshold  is ", threshold);

    return make_unique<simulator_type>(
            delta, threshold, step_lim, save_step,
            read_system<traitsT>(root, 0),
            read_forcefield<traitsT>(root, 0),
            read_observer<traitsT>(root));
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulated_annealing_simulator(
        const toml::table& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type   = typename traitsT::real_type;

    const auto& integrator     = toml::find(simulator, "integrator");
    const auto integrator_type = toml::find<std::string>(integrator, "type");
    const auto tstep      = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep      = toml::find<std::size_t>(simulator, "save_step");

    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    const auto schedule  = toml::find(simulator, "schedule");
    const auto schedule_type  = toml::find<std::string>(schedule, "type");
    const auto schedule_begin = toml::find<real_type>  (schedule, "begin");
    const auto schedule_end   = toml::find<real_type>  (schedule, "end");
    const auto each_step = toml::find<std::size_t>(simulator, "each_step");

    MJOLNIR_LOG_NOTICE("temperature from ", schedule_begin);
    MJOLNIR_LOG_NOTICE("temperature to   ", schedule_end);
    MJOLNIR_LOG_INFO("update temperature for each ", each_step, " steps");

    if(schedule_type == "linear")
    {
        MJOLNIR_LOG_NOTICE("temparing schedule is linear.");
        if(integrator_type == "VelocityVerlet")
        {
            MJOLNIR_LOG_ERROR("Simulated Annealing + NVE Newtonian");
            MJOLNIR_LOG_ERROR("NVE Newtonian doesn't have temperature control.");
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_simulated_annealing_simulator: invalid integrator: ",
                toml::find(integrator, "type"), "here", {
                "Newtonian Integrator does not controls temperature."
                "expected value is one of the following.",
                "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                          " based on the Velocity Verlet"
                }));
        }
        else if(integrator_type == "UnderdampedLangevin")
        {
            using integrator_t = UnderdampedLangevinIntegrator<traitsT>;
            using simulator_t  = SimulatedAnnealingSimulator<
                traitsT, integrator_t, LinearScheduler>;

            MJOLNIR_LOG_NOTICE("Integrator is Underdamped Langevin.");
            return make_unique<simulator_t>(tstep, sstep, each_step,
                    LinearScheduler<real_type>(schedule_begin, schedule_end),
                    read_system<traitsT>(root, 0),
                    read_forcefield<traitsT>(root, 0),
                    read_underdamped_langevin_integrator<traitsT>(simulator),
                    read_observer<traitsT>(root));
        }
        else if(integrator_type == "BAOABLangevin")
        {
            using integrator_t = BAOABLangevinIntegrator<traitsT>;
            using simulator_t  = SimulatedAnnealingSimulator<
                traitsT, integrator_t, LinearScheduler>;

            MJOLNIR_LOG_NOTICE("Integrator is BAOAB Langevin.");
            return make_unique<simulator_t>(tstep, sstep, each_step,
                    LinearScheduler<real_type>(schedule_begin, schedule_end),
                    read_system<traitsT>(root, 0),
                    read_forcefield<traitsT>(root, 0),
                    read_BAOAB_langevin_integrator<traitsT>(simulator),
                    read_observer<traitsT>(root));
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_simulated_annealing_simulator: invalid integrator: ",
                toml::find(simulator, "integrator"), "here", {
                "expected value is one of the following.",
                "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                          " based on the Velocity Verlet"
                }));
        }
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_simulated_annealing_simulator: invalid schedule.type",
            toml::find<toml::value>(schedule, "type"), "here", {
            "expected value is one of the following.",
            "- \"linear\"     : simple linear temperature scheduling",
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulator_from_table(const toml::table& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto type = toml::find<std::string>(simulator, "type");
    if(type == "MolecularDynamics")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is Molecular Dynamics.");
        return read_molecular_dynamics_simulator<traitsT>(root, simulator);
    }
    else if(type == "SteepestDescent")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is Steepest Descent.");
        return read_steepest_descent_simulator<traitsT>(root, simulator);
    }
    else if(type == "SimulatedAnnealing")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is Simulated Annealing.");
        return read_simulated_annealing_simulator<traitsT>(root, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_simulator: invalid type",
            toml::find<toml::value>(simulator, "type"), "here", {
            "expected value is one of the following.",
            "- \"MolecularDynamcis\" : standard MD simulation",
            "- \"SteepestDescent\"   : energy minimization by gradient method",
            "- \"SimulatedAnnealing\": energy minimization by Annealing",
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulator(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto& simulator  = toml::find(root, "simulator");
    if(simulator.as_table().count("file_name") == 1)
    {
        MJOLNIR_LOG_SCOPE(if(simulator.as_table().count("file_name") == 1));

        const auto input_path = read_input_path(root);
        const auto file_name  = toml::find<std::string>(simulator, "file_name");
        MJOLNIR_LOG_INFO("file_name = ", file_name);

        if(simulator.as_table().size() != 1)
        {
            MJOLNIR_LOG_WARN("[simulator] has `file_name` key and other keys.");
            MJOLNIR_LOG_WARN("When `file_name` is provided, other values are "
                             "ignored because those are read from the specified"
                             " file (", input_path, file_name, ").");
        }

        MJOLNIR_LOG_NOTICE("simulator is defined in ", input_path, file_name);
        MJOLNIR_LOG_NOTICE("reading ", input_path, file_name, " ...");
        const auto simfile = toml::parse(input_path + file_name);
        MJOLNIR_LOG_NOTICE(" done.");

        if(simfile.count("simulator") != 1)
        {
            throw_exception<std::out_of_range>("[error] mjolnir::read_simulator: "
                "table [simulator] not found in the toml file\n --> ",
                input_path, file_name, "\n | the file should define [simulator] "
                "table and define values in it.");
        }
        return read_simulator_from_table<traitsT>(root, simfile.at("simulator"));
    }
    else
    {
        return read_simulator_from_table<traitsT>(root, simulator);
    }
}

} // mjolnir
#endif// MJOLNIR_READ_SIMULATOR
