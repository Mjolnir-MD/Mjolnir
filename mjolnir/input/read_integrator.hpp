#ifndef MJOLNIR_READ_INTEGRATOR
#define MJOLNIR_READ_INTEGRATOR
#include <extlib/toml/toml/toml.hpp>
#include <mjolnir/core/VelocityVerletStepper.hpp>
#include <mjolnir/core/UnderdampedLangevinStepper.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
VelocityVerletStepper<traitsT>
read_velocity_verlet_stepper(const toml::table& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_velocity_verlet_stepper(), 0);
    typedef typename traitsT::real_type real_type;

    const real_type delta_t =
        get_toml_value<real_type>(simulator, "delta_t", "[simulator]");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    return VelocityVerletStepper<traitsT>(delta_t);
}


template<typename traitsT>
UnderdampedLangevinStepper<traitsT>
read_underdamped_langevin_stepper(const toml::table& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_underdamped_langevin_stepper(), 0);
    typedef typename traitsT::real_type real_type;

    const auto& parameters = // array of tables
        get_toml_value<toml::array>(simulator, "parameters", "[simulator]");
    const auto  seed =
        get_toml_value<std::uint32_t>(simulator, "seed", "[simulator]");

    std::vector<real_type> gamma(parameters.size());
    for(const auto& tab : parameters)
    {
        const auto& params = toml::get<toml::table>(tab);
        const std::size_t idx = get_toml_value<std::size_t>(
                params, "index", "<anonymous> in register");
        const real_type    gm = get_toml_value<real_type>(
                params, "gamma", "<anonymous> in register");

        if(gamma.size() <= idx){gamma.resize(idx+1);}
        gamma.at(idx) = gm;
        MJOLNIR_LOG_INFO("idx = ", idx, ", gamma = ", gm);
    }

    const real_type delta_t =
        get_toml_value<real_type>(simulator, "delta_t", "[simulator]");
    MJOLNIR_LOG_INFO("delta_t = ", delta_t);

    return UnderdampedLangevinStepper<traitsT>(delta_t, std::move(gamma),
            RandomNumberGenerator<traitsT>(seed));
}

}
#endif// MJOLNIR_READ_INTEGRATOR
