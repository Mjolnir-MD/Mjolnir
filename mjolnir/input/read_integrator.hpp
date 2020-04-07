#ifndef MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#define MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/throw_exception.hpp>

namespace mjolnir
{

template<typename traitsT>
SystemMotionRemover<traitsT>
read_system_motion_remover(const toml::value& simulator)
{
    if(!simulator.contains("integrator") ||
       !simulator.at("integrator").contains("remove"))
    {
        return SystemMotionRemover<traitsT>(false, false, false);
    }

    const auto& remove = toml::find(simulator, "integrator", "remove");

    const bool translation = toml::find_or(remove, "translation", false);
    const bool rotation    = toml::find_or(remove, "rotation",    false);
    const bool rescale     = toml::find_or(remove, "rescale",     false);

    return SystemMotionRemover<traitsT>(translation, rotation, rescale);
}


template<typename traitsT>
VelocityVerletIntegrator<traitsT>
read_velocity_verlet_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    return VelocityVerletIntegrator<traitsT>(delta_t,
            read_system_motion_remover<traitsT>(simulator));
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

    check_keys_available(integrator, {"type"_s, "seed"_s, "parameters"_s, "remove"_s, "env"_s});

    const auto parameters = toml::find<toml::array  >(integrator, "parameters");
    const auto& env = simulator.as_table().count("env") == 1 ?
                      simulator.as_table().at("env") : toml::value{};

    std::vector<real_type> gamma(parameters.size());
    for(const auto& params : parameters)
    {
        const auto offset = find_parameter_or<std::int64_t>(params, env, "offset", 0);
        const auto idx = toml::find<std::size_t>(params, "index") + offset;
        const auto  gm = find_parameter<real_type>(params, env, "gamma", u8"γ");
        if(gamma.size() <= idx){gamma.resize(idx+1);}
        gamma.at(idx) = gm;

        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }
    return UnderdampedLangevinIntegrator<traitsT>(delta_t, std::move(gamma),
            read_system_motion_remover<traitsT>(simulator));
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

    check_keys_available(integrator, {"type"_s, "seed"_s, "parameters"_s, "remove"_s, "env"_s});

    const auto parameters = toml::find<toml::array  >(integrator, "parameters");
    const auto& env = simulator.as_table().count("env") == 1?
                      simulator.as_table().count("env") : toml::value{};

    std::vector<real_type> gamma(parameters.size());
    for(const auto& params : parameters)
    {
        const auto offset = find_parameter_or<std::int64_t>(params, env, "offset", 0);
        const auto idx = find_parameter<std::size_t>(params, env, "index") + offset;
        const auto  gm = find_parameter<real_type>(params, env, "gamma", u8"γ");
        if(gamma.size() <= idx) {gamma.resize(idx+1);}
        gamma.at(idx) = gm;

        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }
    return BAOABLangevinIntegrator<traitsT>(delta_t, std::move(gamma),
            read_system_motion_remover<traitsT>(simulator));
}

// A mapping object from type information (template parameter) to the actual
// read_xxx_integrator function
template<typename T>
struct read_integrator_impl;

template<typename traitsT>
struct read_integrator_impl<VelocityVerletIntegrator<traitsT>>
{
    static VelocityVerletIntegrator<traitsT> invoke(const toml::value& sim)
    {
        return read_velocity_verlet_integrator<traitsT>(sim);
    }
};

template<typename traitsT>
struct read_integrator_impl<UnderdampedLangevinIntegrator<traitsT>>
{
    static UnderdampedLangevinIntegrator<traitsT> invoke(const toml::value& sim)
    {
        return read_underdamped_langevin_integrator<traitsT>(sim);
    }
};

template<typename traitsT>
struct read_integrator_impl<BAOABLangevinIntegrator<traitsT>>
{
    static BAOABLangevinIntegrator<traitsT> invoke(const toml::value& sim)
    {
        return read_BAOAB_langevin_integrator<traitsT>(sim);
    }
};

template<typename integratorT>
integratorT read_integrator(const toml::value& sim)
{
    if(sim.as_table().count("integrator") == 0)
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_integrator: No integrator defined: ", sim, "here", {
            "expected value is one of the following.",
            "- \"VelocityVerlet\"     : simple and standard Velocity Verlet integrator.",
            "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                      " based on the Velocity Verlet",
            "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator"
            }));
    }
    return read_integrator_impl<integratorT>::invoke(sim);
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_velocity_verlet_integrator(const toml::value& simulator);
extern template VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_velocity_verlet_integrator(const toml::value& simulator);
extern template VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_velocity_verlet_integrator(const toml::value& simulator);
extern template VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_velocity_verlet_integrator(const toml::value& simulator);

extern template UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_underdamped_langevin_integrator(const toml::value& simulator);
extern template UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_underdamped_langevin_integrator(const toml::value& simulator);
extern template UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_underdamped_langevin_integrator(const toml::value& simulator);
extern template UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_underdamped_langevin_integrator(const toml::value& simulator);

extern template BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_BAOAB_langevin_integrator(const toml::value& simulator);
extern template BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_BAOAB_langevin_integrator(const toml::value& simulator);
extern template BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_BAOAB_langevin_integrator(const toml::value& simulator);
extern template BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_BAOAB_langevin_integrator(const toml::value& simulator);
#endif

} // mjolnir
#endif// MJOLNIR_READ_INTEGRATOR
