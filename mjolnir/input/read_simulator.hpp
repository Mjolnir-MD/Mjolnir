#ifndef MJOLNIR_READ_SIMULATOR
#define MJOLNIR_READ_SIMULATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/MDSimulator.hpp>
#include <mjolnir/core/SteepestDescentSimulator.hpp>
#include <mjolnir/core/SimulatedAnnealingSimulator.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/get_toml_value.hpp>
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

    const std::string integration = toml::get<std::string>(toml_value_at(
            simulator, "scheme", "[simulator]"));
    const std::size_t tstep = toml::get<std::size_t>(toml_value_at(
            simulator, "total_step", "[simulator]"));
    const std::size_t sstep = toml::get<std::size_t>(toml_value_at(
            simulator, "save_step", "[simulator]"));

    MJOLNIR_LOG_INFO("total step = ", tstep);

    if(integration == "Newtonian")
    {
        MJOLNIR_SCOPE(integration == "Newtonian", 2);
        using integrator_t = VelocityVerletStepper<traitsT>;
        using simulator_t  = MDSimulator<traitsT, integrator_t>;
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
        using simulator_t  = MDSimulator<traitsT, integrator_t>;
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

    const std::size_t step_lim  = toml::get<std::size_t>(toml_value_at(
            simulator, "step_limit", "[simulator]"));
    const std::size_t save_step = toml::get<std::size_t>(toml_value_at(
            simulator, "save_step", "[simulator]"));
    const real_type   delta     = toml::get<real_type>(toml_value_at(
            simulator, "delta", "[simulator]"));
    const real_type   threshold = toml::get<real_type>(toml_value_at(
            simulator, "threshold", "[simulator]"));

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

    const std::string integration = toml::get<std::string>(toml_value_at(
            simulator, "scheme", "[simulator]"));
    const std::size_t tstep = toml::get<std::size_t>(toml_value_at(
            simulator, "total_step", "[simulator]"));
    const std::size_t sstep = toml::get<std::size_t>(toml_value_at(
            simulator, "save_step", "[simulator]"));

    MJOLNIR_LOG_INFO("total step = ", tstep);
    MJOLNIR_LOG_INFO("save  step = ", sstep);

    const std::string schedule = toml::get<std::string>(toml_value_at(
            simulator, "schedule", "[simulator]"));
    const real_type   T_from = toml::get<real_type>(toml_value_at(
            simulator, "T_begin", "[simulator]"));
    const real_type   T_to   = toml::get<real_type>(toml_value_at(
            simulator, "T_end",  "[simulator]"));
    const std::size_t each_step = toml::get<std::size_t>(toml_value_at(
            simulator, "each_step",  "[simulator]"));

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

    const std::string type = toml::get<std::string>(
            toml_value_at(simulator, "type", "[simulator]"));

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
    const auto& simulator = toml_value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();

    if(simulator.count("file_name") == 1)
    {
        MJOLNIR_SCOPE(simulator.count("file_name") == 1, 1);
        if(simulator.size() != 1)
        {
            std::cerr << "WARNING: [simulator] has `file_name` key.\n";
            std::cerr << "       : When `file_name` is provided, all settings ";
            std::cerr << "are read from the file, so other fields are ignored.";
            std::cerr << std::endl;
            MJOLNIR_LOG_WARN("[simulator] has file_name and other settings");
        }

        const std::string file_name =
            toml::get<std::string>(simulator.at("file_name"));
        MJOLNIR_LOG_INFO("file_name = ", file_name);

        const auto simulator_file = toml::parse(file_name);
        if(simulator_file.count("simulator") == 1)
        {
            MJOLNIR_LOG_ERROR("`simulator` value found in ", file_name);

            const auto simulator_toml_type =
                simulator_file.at("simulator").type();
            if(simulator_toml_type != toml::value_t::Table)
            {
                std::cerr << "FATAL: each [simulator] should be provided as ";
                std::cerr << "a table in each file (" << file_name <<  ").\n";
                std::cerr << "       : note: [[...]] means array-of-table. ";
                std::cerr << "please take care.\n";
                std::exit(1);
            }
            std::cerr << "WARNING: in `simulator` file, [simulator] table ";
            std::cerr << "is not necessary.\n";

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
