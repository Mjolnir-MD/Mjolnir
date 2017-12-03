#ifndef MJOLNIR_READ_INTEGRATOR
#define MJOLNIR_READ_INTEGRATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletStepper.hpp>
#include <mjolnir/core/UnderdampedLangevinStepper.hpp>
#include <mjolnir/input/get_toml_value.hpp>

namespace mjolnir
{

template<typename traitsT>
VelocityVerletStepper<traitsT>
read_velocity_verlet_stepper(const toml::Table& data)
{
    typedef typename traitsT::real_type real_type;

    const auto& simulator = detail::value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    return VelocityVerletStepper<traitsT>(toml::get<real_type>(
            detail::value_at(simulator, "delta_t", "[[simulator]]")));
}


template<typename traitsT>
UnderdampedLangevinStepper<traitsT>
read_underdamped_langevin_stepper(const toml::Table& data)
{
    typedef typename traitsT::real_type real_type;

    const auto& simulator   = detail::value_at(data, "simulator", "<root>"
            ).cast<toml::value_t::Table>();
    const auto& registrator = detail::value_at(
            simulator, "register", "[[simulator]]").cast<toml::value_t::Array>();

    const std::uint32_t seed = toml::get<std::uint32_t>(
            detail::value_at(simulator, "seed", "[[simulator]]"));

    std::vector<real_type> gamma(registrator.size());
    for(const auto& tab : registrator)
    {
        const auto& params = tab.cast<toml::value_t::Table>();
        const std::size_t idx = toml::get<std::size_t>(detail::value_at(
                params, "index", "<anonymous> in register"));
        const real_type    gm = toml::get<real_type>(detail::value_at(
                params, "gamma", "<anonymous> in register"));

        if(gamma.size() <= idx){gamma.resize(idx+1);}
        gamma.at(idx) = gm;
    }

    return UnderdampedLangevinStepper<traitsT>(toml::get<real_type>(
            detail::value_at(simulator, "delta_t", "[[simulator]]")),
            std::move(gamma), RandomNumberGenerator<traitsT>(seed));
}

}
#endif// MJOLNIR_READ_INTEGRATOR
