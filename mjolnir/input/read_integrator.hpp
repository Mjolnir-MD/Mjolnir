#ifndef MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#define MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/read_utility.hpp>
#include <mjolnir/core/VelocityVerletIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
VelocityVerletIntegrator<traitsT>
read_velocity_verlet_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type = typename traitsT::real_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    return VelocityVerletIntegrator<traitsT>(delta_t);
}

template<typename traitsT>
UnderdampedLangevinIntegrator<traitsT>
read_underdamped_langevin_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    const auto& integrator = toml::find(simulator, "integrator");

    const auto seed         = toml::find<std::uint32_t>(integrator, "seed");
    const auto gamma_reader = [](const toml::value& v) -> real_type
    {
        const auto& tab = toml::get<toml::table>(v);
        if(tab.count(u8"γ") == 1)
        {
            return toml::find<real_type>(v, u8"γ");
        }
        return toml::find<real_type>(v, "gamma");
    };

    auto parameters = read_array<real_type>(
            toml::find(integrator, "parameters"), gamma_reader);

    return UnderdampedLangevinIntegrator<traitsT>(
            delta_t, std::move(parameters),
            RandomNumberGenerator<traitsT>(seed));
}

} // mjolnir
#endif// MJOLNIR_READ_INTEGRATOR
