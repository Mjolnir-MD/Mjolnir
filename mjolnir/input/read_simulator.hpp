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
read_simulator(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_simulator(), 0);

    using real_type   = typename traitsT::real_type;
    const auto& simulator = toml_value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    const std::string type = toml::get<std::string>(
            toml_value_at(simulator, "type", "[simulator]"));

    if(type == "Molecular Dynamics")
    {
        MJOLNIR_SCOPE(type == "Molecular Dynamics", 1);
        const std::string integration = toml::get<std::string>(toml_value_at(
                simulator, "scheme", "[simulator]"));
        const std::size_t tstep = toml::get<std::size_t>(toml_value_at(
                simulator, "total_step", "[simulator]"));

        MJOLNIR_LOG_INFO("total step = ", tstep);

        if(integration == "Newtonian")
        {
            MJOLNIR_SCOPE(integration == "Newtonian", 2);
            using integrator_t = VelocityVerletStepper<traitsT>;
            using simulator_t  = MDSimulator<traitsT, integrator_t>;
            return make_unique<simulator_t>(
                    tstep,
                    read_system<traitsT>(data, 0),
                    read_forcefield<traitsT>(data, 0),
                    read_velocity_verlet_stepper<traitsT>(data),
                    read_observer<traitsT>(data));
        }
        else if(integration == "Underdamped Langevin")
        {
            MJOLNIR_SCOPE(integration == "Underdamped Langevin", 2);
            using integrator_t = UnderdampedLangevinStepper<traitsT>;
            using simulator_t  = MDSimulator<traitsT, integrator_t>;
            return make_unique<simulator_t>(
                    tstep,
                    read_system<traitsT>(data, 0),
                    read_forcefield<traitsT>(data, 0),
                    read_underdamped_langevin_stepper<traitsT>(data),
                    read_observer<traitsT>(data));
        }
        else
        {
            throw_exception<std::runtime_error>("invalid integration scheme: ",
                    integration, " for simulator ", type);
        }
    }
    else if(type == "Steepest Descent")
    {
        MJOLNIR_SCOPE(type == "Steepest Descent", 1);
        using simulator_t = SteepestDescentSimulator<traitsT>;

        const std::size_t step_lim  = toml::get<std::size_t>(toml_value_at(
                simulator, "step_limit", "[simulator]"));
        const real_type   delta     = toml::get<real_type>(toml_value_at(
                simulator, "delta", "[simulator]"));
        const real_type   threshold = toml::get<real_type>(toml_value_at(
                simulator, "threshold", "[simulator]"));

        MJOLNIR_LOG_INFO("step_lim  = ", step_lim);
        MJOLNIR_LOG_INFO("delta     = ", delta);
        MJOLNIR_LOG_INFO("threshold = ", threshold);

        return make_unique<simulator_t>(
                delta, threshold, step_lim,
                read_system<traitsT>(data, 0),
                read_forcefield<traitsT>(data, 0),
                read_observer<traitsT>(data));
    }
    else if(type == "Simulated Annealing")
    {
        MJOLNIR_SCOPE(type == "Simulated Annealing", 1);
        const std::string integration = toml::get<std::string>(toml_value_at(
                simulator, "scheme", "[simulator]"));
        const std::size_t tstep = toml::get<std::size_t>(toml_value_at(
                simulator, "total_step", "[simulator]"));

        MJOLNIR_LOG_INFO("total step = ", tstep);

        const std::string schedule = toml::get<std::string>(toml_value_at(
                simulator, "schedule", "[simulator]"));
        const real_type   T_from = toml::get<real_type>(toml_value_at(
                simulator, "temperature_from", "[simulator]"));
        const real_type   T_to   = toml::get<real_type>(toml_value_at(
                simulator, "temperature_to",  "[simulator]"));
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
                return make_unique<simulator_t>(tstep, each_step,
                        linear_schedule<real_type>(T_from, T_to),
                        read_system<traitsT>(data, 0),
                        read_forcefield<traitsT>(data, 0),
                        read_underdamped_langevin_stepper<traitsT>(data),
                        read_observer<traitsT>(data));
            }
            else
            {
                throw_exception<std::runtime_error>("invalid integration scheme: ",
                        integration, " for simulator ", type);
            }
        }
    }
    else
    {
        throw_exception<std::runtime_error>("invalid simulator type: ", type);
    }
}

}
#endif// MJOLNIR_READ_SIMULATOR
