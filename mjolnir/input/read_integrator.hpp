#ifndef MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#define MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
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

    const auto seed       = toml::find<std::uint32_t>(integrator, "seed");
    const auto parameters = toml::find<toml::array  >(integrator, "parameters");

    std::vector<real_type> gamma(parameters.size());
    for(const auto& params : parameters)
    {
        const auto idx = toml::find<std::size_t>(params, "index");
        const auto  gm = toml::expect<real_type>(params, u8"γ").or_other(
                         toml::expect<real_type>(params, "gamma")).unwrap();
        if(gamma.size() <= idx){gamma.resize(idx+1);}
        gamma.at(idx) = gm;

        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }
    return UnderdampedLangevinIntegrator<traitsT>(delta_t, std::move(gamma),
            RandomNumberGenerator<traitsT>(seed));
}

template<typename traitsT>
BAOABLangevinIntegrator<traitsT>
read_BAOAB_langevin_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    const auto& integrator = toml::find(simulator, "integrator");

    const auto seed       = toml::find<std::uint32_t>(integrator, "seed");
    const auto parameters = toml::find<toml::array  >(integrator, "parameters");

    std::vector<real_type> gamma(parameters.size());
    for(const auto& params : parameters)
    {
        const auto idx = toml::find<std::size_t>(params, "index");
        const auto  gm = toml::expect<real_type>(params, u8"γ").or_other(
                         toml::expect<real_type>(params, "gamma")).unwrap();
        if(gamma.size() <= idx) {gamma.resize(idx+1);}
        gamma.at(idx) = gm;

        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }
    return BAOABLangevinIntegrator<traitsT>(delta_t, std::move(gamma),
            RandomNumberGenerator<traitsT>(seed));
}

} // mjolnir
#endif// MJOLNIR_READ_INTEGRATOR
