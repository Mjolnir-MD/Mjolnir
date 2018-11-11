#ifndef MJOLNIR_READ_SIMULATOR
#define MJOLNIR_READ_SIMULATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/MolecularDynamicsSimulator.hpp>
#include <mjolnir/core/SteepestDescentSimulator.hpp>
#include <mjolnir/core/SimulatedAnnealingSimulator.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_system.hpp>
#include <mjolnir/input/read_forcefield.hpp>
#include <mjolnir/input/read_integrator.hpp>
#include <mjolnir/input/read_observer.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_molecular_dynamics_simulator(
        const toml::Table& data, const toml::Table& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_molecular_dynamics_simulator(), 0);
    using real_type = typename traitsT::real_type;

    const std::string integration = get_toml_value<std::string>(
            simulator, "scheme", "[simulator]");
    const std::size_t tstep = get_toml_value<std::size_t>(
            simulator, "total_step", "[simulator]");
    const std::size_t sstep = get_toml_value<std::size_t>(
            simulator, "save_step",  "[simulator]");

    MJOLNIR_LOG_INFO("total step = ", tstep);

    if(integration == "Newtonian")
    {
        MJOLNIR_SCOPE(integration == "Newtonian", 2);
        using integrator_t = VelocityVerletStepper<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;
        return make_unique<simulator_t>(
                tstep, sstep,
                read_system<traitsT>(data, 0),
                read_forcefield<traitsT>(data, 0),
                read_velocity_verlet_stepper<traitsT>(simulator),
                read_observer<traitsT>(data));
    }
    else if(integration == "Underdamped Langevin")
    {
        MJOLNIR_SCOPE(integration == "Underdamped Langevin", 2);
        using integrator_t = UnderdampedLangevinStepper<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;
        return make_unique<simulator_t>(
                tstep, sstep,
                read_system<traitsT>(data, 0),
                read_forcefield<traitsT>(data, 0),
                read_underdamped_langevin_stepper<traitsT>(simulator),
                read_observer<traitsT>(data));
    }
    else
    {
        throw_exception<std::runtime_error>("invalid integration scheme: ",
                integration, " for MolecularDynamicsSimulator");
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_steepest_descent_simulator(
        const toml::Table& data, const toml::Table& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_steepest_descent_simulator(), 0);
    using real_type = typename traitsT::real_type;
    using simulator_type = SteepestDescentSimulator<traitsT>;

    const std::size_t step_lim  = get_toml_value<std::size_t>(
            simulator, "step_limit", "[simulator]");
    const std::size_t save_step = get_toml_value<std::size_t>(
            simulator, "save_step", "[simulator]");
    const real_type   delta     = get_toml_value<real_type>(
            simulator, "delta", "[simulator]");
    const real_type   threshold = get_toml_value<real_type>(
            simulator, "threshold", "[simulator]");

    MJOLNIR_LOG_INFO("step_lim  = ", step_lim);
    MJOLNIR_LOG_INFO("save_step = ", save_step);
    MJOLNIR_LOG_INFO("delta     = ", delta);
    MJOLNIR_LOG_INFO("threshold = ", threshold);

    return make_unique<simulator_type>(
            delta, threshold, step_lim, save_step,
            read_system<traitsT>(data, 0),
            read_forcefield<traitsT>(data, 0),
            read_observer<traitsT>(data));
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulated_annealing_simulator(
        const toml::Table& data, const toml::Table& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_simulated_annealing_simulator(), 0);
    using real_type   = typename traitsT::real_type;

    const std::string integration = get_toml_value<std::string>(
            simulator, "scheme", "[simulator]");
    const std::size_t tstep = get_toml_value<std::size_t>(
            simulator, "total_step", "[simulator]");
    const std::size_t sstep = get_toml_value<std::size_t>(
            simulator, "save_step", "[simulator]");

    MJOLNIR_LOG_INFO("total step = ", tstep);
    MJOLNIR_LOG_INFO("save  step = ", sstep);

    const std::string schedule = get_toml_value<std::string>(
            simulator, "schedule", "[simulator]");
    const real_type   T_from = get_toml_value<real_type>(
            simulator, "T_begin", "[simulator]");
    const real_type   T_to   = get_toml_value<real_type>(
            simulator, "T_end",  "[simulator]");
    const std::size_t each_step = get_toml_value<std::size_t>(
            simulator, "each_step",  "[simulator]");

    MJOLNIR_LOG_INFO("temperature from = ", T_from);
    MJOLNIR_LOG_INFO("temperature to   = ", T_to);
    MJOLNIR_LOG_INFO("for each step    = ", each_step);

    if(schedule == "linear")
    {
        if(integration == "Newtonian")
        {
            MJOLNIR_LOG_ERROR("Simulated Annealing + NVE Newtonian");
            MJOLNIR_LOG_ERROR("NVE Newtonian doesn't have temperature control.");
            throw_exception<std::runtime_error>("Simulated Annealing has ",
                    "no effect for Newtonian Integrator");
        }
        else if(integration == "Underdamped Langevin")
        {
            using integrator_t = UnderdampedLangevinStepper<traitsT>;
            using simulator_t  = SimulatedAnnealingSimulator<
                traitsT, integrator_t, linear_schedule>;

            MJOLNIR_LOG_INFO("Underdamped Langevin is used as a MD engine");
            return make_unique<simulator_t>(tstep, sstep, each_step,
                    linear_schedule<real_type>(T_from, T_to),
                    read_system<traitsT>(data, 0),
                    read_forcefield<traitsT>(data, 0),
                    read_underdamped_langevin_stepper<traitsT>(simulator),
                    read_observer<traitsT>(data));
        }
        else
        {
            throw_exception<std::runtime_error>("invalid integration scheme: ",
                    integration, " for Simulated Annealing");
        }
    }
    else
    {
        throw_exception<std::runtime_error>("invalid temperature schedule: ",
                schedule, " for Simulated Annealing");
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulator_from_table(const toml::Table& data, const toml::Table& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_simulator_from_table(), 0);

    const std::string type =
        get_toml_value<std::string>(simulator, "type", "[simulator]");

    if(type == "Molecular Dynamics")
    {
        MJOLNIR_SCOPE(type == "Molecular Dynamics", 1);
        return read_molecular_dynamics_simulator<traitsT>(data, simulator);
    }
    else if(type == "Steepest Descent")
    {
        MJOLNIR_SCOPE(type == "Steepest Descent", 1);
        return read_steepest_descent_simulator<traitsT>(data, simulator);
    }
    else if(type == "Simulated Annealing")
    {
        MJOLNIR_SCOPE(type == "Simulated Annealing", 1);
        return read_simulated_annealing_simulator<traitsT>(data, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>("invalid simulator type: ", type);
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulator(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_simulator(), 0);

    using real_type   = typename traitsT::real_type;
    const auto& simulator =
        get_toml_value<toml::Table>(data, "simulator", "<root>");

    const auto& files = get_toml_value<toml::Table>(data, "files", "<root>");
    std::string input_path_("./");
    if(files.count("input_path") == 1)
    {
        input_path_ = get_toml_value<std::string>(files, "input_path", "[files]");
    }
    const auto input_path(input_path_); // add constness

    if(simulator.count("file_name") == 1)
    {
        MJOLNIR_SCOPE(simulator.count("file_name") == 1, 1);

        const std::string file_name =
            get_toml_value<std::string>(simulator, "file_name", "[simulator]");
        MJOLNIR_LOG_INFO("file_name = ", file_name);

        if(simulator.size() != 1)
        {
            MJOLNIR_LOG_WARN("[simulator] has `file_name` key and other keys.");
            MJOLNIR_LOG_WARN("When `file_name` is provided, other values are "
                             "ignored because those are read from the specified"
                             " file (", file_name, ").");
        }

        const auto simulator_file = toml::parse(input_path + file_name);
        if(simulator_file.count("simulator") == 1)
        {
            MJOLNIR_LOG_WARN("in `simulator` file, root object is treated as "
                             "a [simulator] table.");
            MJOLNIR_LOG_WARN("but in ", file_name, ", `simulator` key found. "
                             "trying to read it as a simulator setup.");

            if(simulator_file.at("simulator").type() != toml::value_t::Table)
            {
                MJOLNIR_LOG_ERROR("type of `simulator` is different from "
                                  "toml::Table in file (", file_name, ").");
                MJOLNIR_LOG_ERROR("note: [[...]] means Array-of-Tables. "
                                  "please take care.");
                std::exit(1);
            }

            MJOLNIR_LOG_INFO("reading `[simulator]` table");
            return read_simulator_from_table<traitsT>(data, simulator_file.at(
                    "simulator").template cast<toml::value_t::Table>());
        }
        return read_simulator_from_table<traitsT>(data, simulator_file);
    }
    else
    {
        return read_simulator_from_table<traitsT>(data, simulator);
    }
}

}
#endif// MJOLNIR_READ_SIMULATOR
