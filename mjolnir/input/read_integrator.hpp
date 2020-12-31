#ifndef MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#define MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
#include <mjolnir/core/gBAOABLangevinIntegrator.hpp>
#include <mjolnir/core/GJFNVTLangevinIntegrator.hpp>
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

    check_keys_available(integrator, {"type"_s, "seed"_s, "remove"_s, "env"_s,
                                      "parameters"_s, "gammas"_s});

    toml::array gammas;
    if(integrator.contains("parameters"))
    {
        if(integrator.contains("gammas"))
        {
            throw std::runtime_error(toml::format_error("duplicated keys: both "
                "\"parameters\" and \"gammas\" are provided"_s,
                integrator.at("gammas"), "\"gammas\" is defined here"_s,
                integrator.at("parameters"), "\"parameters\" is defined here"_s));
        }
        MJOLNIR_LOG_WARN("deprecated: key \"parameters\" in an integrator is "
                         "depricated. use \"gammas\" instead.");
        gammas = toml::find<toml::array>(integrator, "parameters");
    }
    else
    {
        gammas = toml::find<toml::array>(integrator, "gammas");
    }

    const auto& env = integrator.contains("env") ?
                      integrator.at("env") : toml::value{};

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

    check_keys_available(integrator, {"type"_s, "seed"_s, "remove"_s, "env"_s,
                                      "parameters"_s, "gammas"_s});

    toml::array gammas;
    if(integrator.contains("parameters"))
    {
        if(integrator.contains("gammas"))
        {
            throw std::runtime_error(toml::format_error("duplicated keys: both "
                "\"parameters\" and \"gammas\" are provided"_s,
                integrator.at("gammas"), "\"gammas\" is defined here"_s,
                integrator.at("parameters"), "\"parameters\" is defined here"_s));
        }
        MJOLNIR_LOG_WARN("deprecated: key \"parameters\" in an integrator is "
                         "depricated. use \"gammas\" instead.");
        gammas = toml::find<toml::array>(integrator, "parameters");
    }
    else
    {
        gammas = toml::find<toml::array>(integrator, "gammas");
    }

    const auto& env = integrator.contains("env") ?
                      integrator.at("env") : toml::value{};

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
    return BAOABLangevinIntegrator<traitsT>(delta_t, std::move(gamma),
            read_system_motion_remover<traitsT>(simulator));
}

template<typename traitsT>
gBAOABLangevinIntegrator<traitsT>
read_gBAOAB_langevin_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    const auto& integrator = toml::find(simulator, "integrator");
    check_keys_available(integrator,
            {"type"_s, "seed"_s, "gammas"_s, "remove"_s, "env"_s});

    const auto& env = integrator.contains("env") ?
                      integrator.at("env") : toml::value{};

    const auto gammas = toml::find<toml::array>(integrator, "gammas");

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

    return gBAOABLangevinIntegrator<traitsT>(delta_t, std::move(gamma),
            read_system_motion_remover<traitsT>(simulator));
}

template<typename traitsT>
GJFNVTLangevinIntegrator<traitsT>
read_GJFNVT_langevin_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    const auto& integrator = toml::find(simulator, "integrator");

    check_keys_available(integrator, {"type"_s, "seed"_s, "remove"_s, "env"_s,
                                      "parameters"_s, "alphas"_s});

    toml::array alphas;
    if(integrator.contains("parameters"))
    {
        if(integrator.contains("alphas"))
        {
            throw std::runtime_error(toml::format_error("duplicated keys: both "
                "\"parameters\" and \"alphas\" are provided"_s,
                integrator.at("alphas"), "\"alphas\" is defined here"_s,
                integrator.at("parameters"), "\"parameters\" is defined here"_s));
        }
        MJOLNIR_LOG_WARN("deprecated: key \"parameters\" in an integrator is "
                         "depricated. use \"alphas\" instead.");
        alphas = toml::find<toml::array>(integrator, "parameters");
    }
    else
    {
        alphas = toml::find<toml::array>(integrator, "alphas");
    }

    const auto& env = integrator.contains("env") ?
                      integrator.at("env") : toml::value{};

    // Temporarily make a vector of idx alpha pair to check index duplication.
    std::vector<std::pair<std::size_t, real_type>> idx_alphas;
    idx_alphas.reserve(alphas.size());
    for(const auto& params : alphas)
    {
        const auto offset = find_parameter_or<std::int64_t>(params, env, "offset", 0);
        const auto idx    = toml::find<std::size_t>(params, "index") + offset;
        const auto al     = find_parameter<real_type>(params, env, "alpha");

        idx_alphas.emplace_back(idx, al);
    }
    check_parameter_overlap(env, alphas, idx_alphas);

    std::vector<real_type> alpha(alphas.size());
    for(const auto& idx_alpha : idx_alphas)
    {
        const auto idx = idx_alpha.first;
        const auto  al = idx_alpha.second;
        if(alpha.size() <= idx)
        {
            alpha.resize(idx+1, real_type(0.0));
        }
        alpha.at(idx) = al;

        MJOLNIR_LOG_INFO("idx = ", idx, ", alpha = ", al);
    }
    return GJFNVTLangevinIntegrator<traitsT>(delta_t, std::move(alpha),
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
struct read_integrator_impl<gBAOABLangevinIntegrator<traitsT>>
{
    static gBAOABLangevinIntegrator<traitsT> invoke(const toml::value& sim)
    {
        return read_gBAOAB_langevin_integrator<traitsT>(sim);
    }
};

template<typename traitsT>
struct read_integrator_impl<GJFNVTLangevinIntegrator<traitsT>>
{
    static GJFNVTLangevinIntegrator<traitsT> invoke(const toml::value& sim)
    {
        return read_GJFNVT_langevin_integrator<traitsT>(sim);
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
            "- \"g-BAOABLangevin\"    : geodesic BAOAB Langevin Integrator",
            "- \"G-JFLangevin\"       : Verlet-type Langevin Integrator by G-J&F"
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

extern template gBAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_gBAOAB_langevin_integrator(const toml::value& simulator);
extern template gBAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_gBAOAB_langevin_integrator(const toml::value& simulator);
extern template gBAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_gBAOAB_langevin_integrator(const toml::value& simulator);
extern template gBAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_gBAOAB_langevin_integrator(const toml::value& simulator);

extern template GJFNVTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_GJFNVT_langevin_integrator(const toml::value& simulator);
extern template GJFNVTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_GJFNVT_langevin_integrator(const toml::value& simulator);
extern template GJFNVTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_GJFNVT_langevin_integrator(const toml::value& simulator);
extern template GJFNVTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_GJFNVT_langevin_integrator(const toml::value& simulator);
#endif

} // mjolnir
#endif// MJOLNIR_READ_INTEGRATOR
