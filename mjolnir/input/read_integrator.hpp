#ifndef MJOLNIR_READ_INTEGRATOR
#define MJOLNIR_READ_INTEGRATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/NewtonianStepper.hpp>
#include <mjolnir/core/UnderdampedLangevinStepper.hpp>

namespace mjolnir
{

template<typename traitsT>
VelocityVerletStepper<traitsT>
read_velocity_verlet_stepper(const toml::Table& data)
{
    const auto& simulator = data.at("simulator").cast<toml::value_t::Table>();
    return VelocityVerletStepper<traitsT>(
            toml::get<typename traitsT::real_type>(simulator.at("delta_t")));
}


template<typename traitsT>
UnderdampedLangevinStepper<traitsT>
read_underdamped_langevin_stepper(const toml::Table& data)
{
    typedef typename traitsT::real_type real_type;
    const auto& simulator = data.at("simulator").cast<toml::value_t::Table>();
    const auto& registrator =
        simulator.at("register").cast<toml::value_t::Array>();

    const std::uint32_t seed = toml::get<std::uint32_t>(simulator.at("seed"));

    std::vector<real_type> gamma(registrator.size());
    for(const auto& tab : registrator)
    {
        const auto& params = tab.cast<toml::value_t::Table>();
        const std::size_t idx = toml::get<std::size_t>(params.at("index"));
        const real_type    gm = toml::get<real_type>(params.at("gamma"));

        if(gamma.size() <= idx) gamma.resize(idx+1);
        gamma.at(idx) = gm;
    }

    return UnderdampedLangevinStepper<traitsT>(
            toml::get<typename traitsT::real_type>(simulator.at("delta_t")),
            std::move(gamma), RandomNumberGenerator<traitsT>(seed));
}

}
#endif// MJOLNIR_READ_INTEGRATOR
