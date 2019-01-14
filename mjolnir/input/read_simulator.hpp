#ifndef MJOLNIR_READ_SIMULATOR
#define MJOLNIR_READ_SIMULATOR
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
#include <mjolnir/input/read_files_table.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_molecular_dynamics_simulator(
        const toml::table& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_molecular_dynamics_simulator(), 0);
    using real_type = typename traitsT::real_type;

    const auto integrator = toml::find<std::string>(simulator, "integrator");
    const auto tstep      = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep      = toml::find<std::size_t>(simulator, "save_step");
    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    if(integrator == "Newtonian")
    {
        MJOLNIR_LOG_NOTICE("Integrator is Newtonian.");
        using integrator_t = VelocityVerletIntegrator<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;

        return make_unique<simulator_t>(
                tstep, sstep,
                read_system<traitsT>(root, 0),
                read_forcefield<traitsT>(root, 0),
                read_velocity_verlet_integrator<traitsT>(simulator),
                read_observer<traitsT>(root));
    }
    else if(integrator == "Underdamped Langevin")
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
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_molecular_dynamics_simulator: invalid integrator: ",
            toml::find(simulator, "integrator"), "here", {
            "expected value is one of the following.",
            "- \"Newtonian\"           : simple and standard Velocity Verlet integrator.",
            "- \"Underdamped Langevin\": simple Underdamped Langevin Integrator"
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
    MJOLNIR_SCOPE(read_steepest_descent_simulator(), 0);
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
    MJOLNIR_SCOPE(read_simulated_annealing_simulator(), 0);
    using real_type   = typename traitsT::real_type;

    const auto integrator = toml::find<std::string>(simulator, "integrator");
    const auto tstep      = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep      = toml::find<std::size_t>(simulator, "save_step");

    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    const auto schedule  = toml::find<std::string>(simulator, "schedule");
    const auto T_begin   = toml::find<real_type>  (simulator, "T_begin");
    const auto T_end     = toml::find<real_type>  (simulator, "T_end");
    const auto each_step = toml::find<std::size_t>(simulator, "each_step");

    MJOLNIR_LOG_NOTICE("temperature from ", T_begin);
    MJOLNIR_LOG_NOTICE("temperature to   ", T_end);
    MJOLNIR_LOG_INFO("update temperature for each ", each_step, " steps");

    if(schedule == "linear")
    {
        MJOLNIR_LOG_NOTICE("temparing schedule is linear.");
        if(integrator == "Newtonian")
        {
            MJOLNIR_LOG_ERROR("Simulated Annealing + NVE Newtonian");
            MJOLNIR_LOG_ERROR("NVE Newtonian doesn't have temperature control.");
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_simulated_annealing_simulator: invalid integrator: ",
                toml::find(simulator, "integrator"), "here", {
                "Newtonian Integrator does not controls temperature."
                "expected value is one of the following.",
                "- \"Underdamped Langevin\": simple Underdamped Langevin Integrator"
                                           " based on the Velocity Verlet"
                }));
        }
        else if(integrator == "Underdamped Langevin")
        {
            using integrator_t = UnderdampedLangevinIntegrator<traitsT>;
            using simulator_t  = SimulatedAnnealingSimulator<
                traitsT, integrator_t, LinearScheduler>;

            MJOLNIR_LOG_NOTICE("Integrator is Underdamped Langevin.");
            return make_unique<simulator_t>(tstep, sstep, each_step,
                    LinearScheduler<real_type>(T_begin, T_end),
                    read_system<traitsT>(root, 0),
                    read_forcefield<traitsT>(root, 0),
                    read_underdamped_langevin_integrator<traitsT>(simulator),
                    read_observer<traitsT>(root));
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_simulated_annealing_simulator: invalid integrator: ",
                toml::find(simulator, "integrator"), "here", {
                "expected value is one of the following.",
                "- \"Underdamped Langevin\": simple Underdamped Langevin Integrator"
                                           " based on the Velocity Verlet"
                }));
        }
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_simulated_annealing_simulator: invalid schedule",
            toml::find<toml::value>(simulator, "schedule"), "here", {
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
    MJOLNIR_SCOPE(read_simulator_from_table(), 0);

    const auto type = toml::find<std::string>(simulator, "type");
    if(type == "Molecular Dynamics")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is Molecular Dynamics.");
        return read_molecular_dynamics_simulator<traitsT>(root, simulator);
    }
    else if(type == "Steepest Descent")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is Steepest Descent.");
        return read_steepest_descent_simulator<traitsT>(root, simulator);
    }
    else if(type == "Simulated Annealing")
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
            "- \"Molecular Dynamcis\" : standard MD simulation",
            "- \"Steepest Descent\"   : energy minimization by gradient method",
            "- \"Simulated Annealing\": energy minimization by Annealing",
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulator(const toml::table& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_simulator(), 0);
    using real_type = typename traitsT::real_type;

    const auto& simulator  = toml::find(root, "simulator");
    if(toml::get<toml::table>(simulator).count("file_name") == 1)
    {
        MJOLNIR_SCOPE(simulator.count("file_name") == 1, 1);

        const auto input_path = read_input_path(root);
        const auto file_name  = toml::find<std::string>(simulator, "file_name");
        MJOLNIR_LOG_INFO("file_name = ", file_name);

        if(toml::get<toml::table>(simulator).size() != 1)
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
