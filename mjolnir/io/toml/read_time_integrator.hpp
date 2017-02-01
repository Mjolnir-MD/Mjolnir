#ifndef MJOLNIR_TOML_READ_INTEGRATOR
#define MJOLNIR_TOML_READ_INTEGRATOR
#include <mjolnir/core/VelocityVerlet.hpp>
#include <mjolnir/core/UnderdampedLangevin.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<Integrator<traitsT>>
read_time_integrator(const toml::Table& sim,
        const std::shared_ptr<RandomNumberGenerator<traitsT>>& rng)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_time_integrator CALLED");

    const std::string integral =
        toml::get<toml::String>(sim.at("time_integration"));
    MJOLNIR_LOG_INFO("time integration", integral);

    if(integral == "Underdamped Langevin")
    {
        const typename traitsT::real_type delta_t =
            toml::get<toml::Float>(sim.at("delta_t"));
        MJOLNIR_LOG_INFO("delta_t", delta_t);

        const typename traitsT::real_type temperature =
            toml::get<toml::Float>(sim.at("temperature"));
        MJOLNIR_LOG_INFO("temperatrue", temperature);

        const typename traitsT::real_type kB =
            toml::get<toml::Float>(sim.at("kB"));
        MJOLNIR_LOG_INFO("kB", kB);

        const std::size_t num_particles =
            toml::get<toml::Integer>(sim.at("number_of_particle"));
        MJOLNIR_LOG_INFO("num_particles", num_particles);

        std::vector<typename traitsT::real_type> fric =
            toml::get<toml::Array<toml::Float>>(sim.at("friction_constant"));
        MJOLNIR_LOG_INFO("fric.size()", fric.size());

        if(num_particles != fric.size())
            throw std::runtime_error("friction constant not enough");

        MJOLNIR_LOG_DEBUG("read_time_integrator RETURNED");
        return make_unique<UnderdampedLangevin<traitsT>>(
            delta_t, num_particles, temperature, kB, std::move(fric), rng);
    }
    else if(integral == "NVE Newtonian")
    {
        const typename traitsT::real_type delta_t =
            toml::get<toml::Float>(sim.at("delta_t"));
        MJOLNIR_LOG_INFO("delta_t", delta_t);

        const std::size_t num_particles =
            toml::get<toml::Integer>(sim.at("number_of_particle"));
        MJOLNIR_LOG_INFO("num_particles", num_particles);

        MJOLNIR_LOG_DEBUG("read_time_integrator RETURNED");
        return make_unique<VelocityVerlet<traitsT>>(delta_t, num_particles);
    }
    else
        throw std::invalid_argument("unknown time integral method: " + integral);
}


} // mjolnir
#endif /* MJOLNIR_TOML_READ_INTEGRATOR */
