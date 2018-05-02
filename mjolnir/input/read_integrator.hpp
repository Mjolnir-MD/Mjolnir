#ifndef MJOLNIR_READ_INTEGRATOR
#define MJOLNIR_READ_INTEGRATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletStepper.hpp>
#include <mjolnir/core/UnderdampedLangevinStepper.hpp>
#include <mjolnir/util/get_toml_value.hpp>

namespace mjolnir
{

template<typename traitsT>
VelocityVerletStepper<traitsT>
read_velocity_verlet_stepper(const toml::Table& data)
{
    typedef typename traitsT::real_type real_type;

    const auto& simulator = toml_value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    return VelocityVerletStepper<traitsT>(toml::get<real_type>(
            toml_value_at(simulator, "delta_t", "[[simulator]]")));
}


template<typename traitsT>
UnderdampedLangevinStepper<traitsT>
read_underdamped_langevin_stepper(const toml::Table& data)
{
    typedef typename traitsT::real_type real_type;

    const auto& simulator  = toml_value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    const auto& parameters = toml_value_at(
            simulator, "parameters", "[[simulator]]").cast<toml::value_t::Array>();

    const std::uint32_t seed = toml::get<std::uint32_t>(
            toml_value_at(simulator, "seed", "[[simulator]]"));

    std::vector<real_type> gamma(parameters.size());
    for(const auto& tab : parameters)
    {
        const auto& params = tab.cast<toml::value_t::Table>();
        const std::size_t idx = toml::get<std::size_t>(toml_value_at(
                params, "index", "<anonymous> in register"));
        const real_type    gm = toml::get<real_type>(toml_value_at(
                params, "gamma", "<anonymous> in register"));

        if(gamma.size() <= idx){gamma.resize(idx+1);}
        gamma.at(idx) = gm;
    }

    return UnderdampedLangevinStepper<traitsT>(toml::get<real_type>(
            toml_value_at(simulator, "delta_t", "[[simulator]]")),
            std::move(gamma), RandomNumberGenerator<traitsT>(seed));
}

}
#endif// MJOLNIR_READ_INTEGRATOR
