#ifndef MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#define MJOLNIR_INPUT_READ_INTEGRATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
#include <mjolnir/core/GFWNpTLangevinIntegrator.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/throw_exception.hpp>

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

    check_keys_available(integrator, {"type"_s, "parameters"_s});

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
    return UnderdampedLangevinIntegrator<traitsT>(delta_t, std::move(gamma));
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

    check_keys_available(integrator, {"type"_s, "parameters"_s});

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
    return BAOABLangevinIntegrator<traitsT>(delta_t, std::move(gamma));
}

template<typename traitsT>
GFWNpTLangevinIntegrator<traitsT>
read_GFW_NpT_langevin_integrator(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const real_type delta_t = toml::find<real_type>(simulator, "delta_t");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    const auto& integrator = toml::find(simulator, "integrator");

    check_keys_available(integrator, {"type"_s, "chi"_s, "cell_mass"_s,
            "cell_initial_velocity"_s, "cell_gamma"_s, "parameters"_s});

    const auto chi        = toml::find<real_type>(integrator, "chi");
    const auto m_cell     = toml::find<std::array<real_type, 3>>(integrator, "cell_mass");
    const auto gamma_cell = toml::find<std::array<real_type, 3>>(integrator, "cell_gamma");
    const auto v_cell_ini = toml::find_or<std::array<real_type, 3>>(integrator,
        "cell_initial_velocity", std::array<real_type, 3>{{0, 0, 0}});

    const auto& parameters = toml::find<toml::array>(integrator, "parameters");

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
    return GFWNpTLangevinIntegrator<traitsT>(delta_t, chi,
        math::make_coordinate<coordinate_type>(m_cell[0],     m_cell[1],     m_cell[2]),
        math::make_coordinate<coordinate_type>(gamma_cell[0], gamma_cell[1], gamma_cell[2]),
        math::make_coordinate<coordinate_type>(v_cell_ini[0], v_cell_ini[1], v_cell_ini[2]),
        std::move(gamma));
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
struct read_integrator_impl<GFWNpTLangevinIntegrator<traitsT>>
{
    static GFWNpTLangevinIntegrator<traitsT> invoke(const toml::value& sim)
    {
        return read_GFW_NpT_langevin_integrator<traitsT>(sim);
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
            "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator",
            "- \"GFWNpTLangevin\"     : BAOAB based NpT Langevin Integrator"
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

extern template GFWNpTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_GFW_NpT_langevin_integrator(const toml::value& simulator);
extern template GFWNpTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_GFW_NpT_langevin_integrator(const toml::value& simulator);
extern template GFWNpTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_GFW_NpT_langevin_integrator(const toml::value& simulator);
extern template GFWNpTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_GFW_NpT_langevin_integrator(const toml::value& simulator);
#endif

} // mjolnir
#endif// MJOLNIR_READ_INTEGRATOR
