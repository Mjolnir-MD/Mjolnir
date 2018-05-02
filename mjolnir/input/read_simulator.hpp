#ifndef MJOLNIR_READ_SIMULATOR
#define MJOLNIR_READ_SIMULATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/MDSimulator.hpp>
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
    const auto& simulator = toml_value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    const std::string type = toml::get<std::string>(
            toml_value_at(simulator, "type", "[simulator]"));

    if(type == "Molecular Dynamics")
    {
        const std::string integration = toml::get<std::string>(toml_value_at(
                simulator, "scheme", "[simulator]"));
        const std::size_t tstep = toml::get<std::size_t>(toml_value_at(
                simulator, "total_step", "[simulator]"));

        if(integration == "Newtonian")
        {
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
    else
    {
        throw_exception<std::runtime_error>("invalid simulator type: ", type);
    }
}

}
#endif// MJOLNIR_READ_SIMULATOR
