#ifndef MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#define MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
#include <mjolnir/core/GBAOABLangevinIntegrator.hpp>
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

    check_keys_available(integrator,
            {"type"_s, "seed"_s, "parameters"_s, "remove"_s, "env"_s});

    const auto parameters = toml::find<toml::array>(integrator, "parameters");

    const auto& env = integrator.contains("env") ?
                      integrator.at("env") : toml::value{};

    // Temporarily make a vector of idx gamma pair to check index duplication.
    std::vector<std::pair<std::size_t, real_type>> idx_gammas;
    idx_gammas.reserve(parameters.size());
    for(const auto& params : parameters)
    {
        const auto offset = find_parameter_or<std::int64_t>(params, env, "offset", 0);
        const auto idx    = toml::find<std::size_t>(params, "index") + offset;
        const auto gm     = find_parameter<real_type>(params, env, "gamma", u8"γ");

        idx_gammas.emplace_back(idx, gm);
    }
    check_parameter_overlap(env, parameters, idx_gammas);

    std::vector<real_type> gamma(parameters.size());
    for(const auto& idx_gamma : idx_gammas)
    {
        const auto idx = idx_gamma.first;
        const auto gm  = idx_gamma.second;
        if(gamma.size() <= idx)
        {
            gamma.resize(idx+1, real_type(0.0));
        }
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

    check_keys_available(integrator,
            {"type"_s, "seed"_s, "parameters"_s, "remove"_s, "env"_s});

    const auto parameters = toml::find<toml::array>(integrator, "parameters");

    const auto& env = integrator.contains("env") ?
                      integrator.at("env") : toml::value{};

    // Temporarily make a vector of idx gamma pair to check index duplication.
    std::vector<std::pair<std::size_t, real_type>> idx_gammas;
    idx_gammas.reserve(parameters.size());
    for(const auto& params : parameters)
    {
        const auto offset = find_parameter_or<std::int64_t>(params, env, "offset", 0);
        const auto idx    = toml::find<std::size_t>(params, "index") + offset;
        const auto gm     = find_parameter<real_type>(params, env, "gamma", u8"γ");

        idx_gammas.emplace_back(idx, gm);
    }
    check_parameter_overlap(env, parameters, idx_gammas);

    std::vector<real_type> gamma(parameters.size());
    for(const auto& idx_gamma : idx_gammas)
    {
        const auto idx = idx_gamma.first;
        const auto  gm = idx_gamma.second;
        if(gamma.size() <= idx)
        {
            gamma.resize(idx+1, real_type(0.0));
        }
        gamma.at(idx) = gm;

        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }
    return BAOABLangevinIntegrator<traitsT>(delta_t, std::move(gamma),
            read_system_motion_remover<traitsT>(simulator));
}

template<typename traitsT>
GBAOABLangevinIntegrator<traitsT>
read_gBAOAB_langevin_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type                = typename traitsT::real_type;
    using indices_t                = std::array<std::size_t, 2>;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    const auto& integrator = toml::find(simulator, "integrator");
    check_keys_available(integrator,
            {"type"_s, "seed"_s, "parameters"_s, "remove"_s, "env"_s});

    const auto& env = integrator.contains("env") ?
                      integrator.at("env") : toml::value{};

    const auto& parameters = toml::find(integrator, "parameters");
    check_keys_available(parameters, {"gammas"_s, "constraints"_s});

    const auto gammas      = toml::find<toml::array>(parameters, "gammas");
    const auto constraints = toml::find<toml::array>(parameters, "constraints");

    // Read idx gamma pair part in parameters
    // Temporarily make a vector of idx gamma pair to check index duplication.
    std::vector<std::pair<std::size_t, real_type>> idx_gammas;
    idx_gammas.reserve(gammas.size());
    for(const auto& params : gammas)
    {
        const auto offset = find_parameter_or<std::int64_t>(params, env, "offset", 0);
        const auto idx    = toml::find<std::size_t>(params, "index") + offset;
        const auto gm     = find_parameter<real_type>(params, env, "gamma", u8"γ");

        idx_gammas.emplace_back(idx, gm);
    }
    check_parameter_overlap(env, gammas, idx_gammas);

    std::vector<real_type> gamma(gammas.size());
    for(const auto& idx_gamma : idx_gammas)
    {
        const auto idx = idx_gamma.first;
        const auto  gm = idx_gamma.second;
        if(gamma.size() <= idx)
        {
            gamma.resize(idx+1, real_type(0.0));
        }
        gamma.at(idx) = gm;

        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }

    std::vector<std::pair<indices_t, real_type>> indices_v0s;
    indices_v0s.reserve(constraints.size());
    for(const auto& params : constraints)
    {
        const auto offset  = find_parameter_or<std::int64_t>(params, env, "offset", 0);
        auto indices = find_parameter<indices_t>(params, env, "indices");
        for(auto& i : indices) { i += offset;}
        const auto v0      = find_parameter<real_type>(params, env, "v0");

        MJOLNIR_LOG_INFO_NO_LF("idxs = ", indices, ", ");
        indices_v0s.emplace_back(indices, v0);

    }
    // TODO : check parameter overlap in indices_v0s.

    return GBAOABLangevinIntegrator<traitsT>(delta_t, std::move(gamma), std::move(indices_v0s),
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

template<typename traitsT>
struct read_integrator_impl<GBAOABLangevinIntegrator<traitsT>>
{
    static GBAOABLangevinIntegrator<traitsT> invoke(const toml::value& sim)
    {
        return read_gBAOAB_langevin_integrator<traitsT>(sim);
    }
};

template<typename integratorT>
integratorT read_integrator(const toml::value& sim)
{
    if(!sim.contains("integrator"))
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_integrator: No integrator defined: ", sim, "here", {
            "expected value is one of the following.",
            "- \"VelocityVerlet\"     : simple and standard Velocity Verlet integrator.",
            "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                      " based on the Velocity Verlet",
            "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator",
            "- \"g-BAOABLangevin\"     : geodesic BAOAB Langevin Integrator"
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
