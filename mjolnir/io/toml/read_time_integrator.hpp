#ifndef MJOLNIR_TOML_READ_INTEGRATOR
#define MJOLNIR_TOML_READ_INTEGRATOR
#include <mjolnir/core/NVENewtonian.hpp>
#include <mjolnir/core/UnderdampedLangevin.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<Integrator<traitsT>>
read_underdamped_langevin(const toml::Table& sim, const std::string boundary,
        const std::shared_ptr<RandomNumberGenerator<traitsT>>& rng)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_time_integrator CALLED");

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

    if(boundary == "Unlimited")
        return make_unique<UnderdampedLangevin<traitsT>>(
            delta_t, num_particles, temperature, kB, std::move(fric), rng);
    else if(boundary == "Periodic")
        return make_unique<
            UnderdampedLangevin<traitsT, PeriodicBoundaryXYZ<traitsT>>>(
                delta_t, num_particles, temperature, kB, std::move(fric), rng);
    else
        throw std::runtime_error("unknown boundary condition: " + boundary);
}

template<typename traitsT>
std::unique_ptr<Integrator<traitsT>>
read_nve_newtonian(const toml::Table& sim, const std::string boundary)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_time_integrator CALLED");

    const typename traitsT::real_type delta_t =
        toml::get<toml::Float>(sim.at("delta_t"));
    MJOLNIR_LOG_INFO("delta_t", delta_t);

    const std::size_t num_particles =
        toml::get<toml::Integer>(sim.at("number_of_particle"));
    MJOLNIR_LOG_INFO("num_particles", num_particles);

    MJOLNIR_LOG_DEBUG("read_time_integrator RETURNED");

    if(boundary == "Unlimited")
        return make_unique<NVENewtonian<traitsT>>(delta_t, num_particles);
    else if(boundary == "Periodic")
        return make_unique<
            NVENewtonian<traitsT, PeriodicBoundaryXYZ<traitsT>>>(
                delta_t, num_particles);
    else
        throw std::runtime_error("unknown boundary condition: " + boundary);
}


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

    const std::string boundary =
        toml::get<toml::String>(sim.at("boundary"));
    MJOLNIR_LOG_INFO("boundary condition", boundary);

    if(boundary == "Periodic")
    {
        const auto lower = toml::get<toml::Array<toml::Float>>(sim.at("lower"));
        const auto upper = toml::get<toml::Array<toml::Float>>(sim.at("upper"));
        const typename traitsT::coordinate_type low(lower.at(0), lower.at(1), lower.at(2));
        const typename traitsT::coordinate_type upp(upper.at(0), upper.at(1), upper.at(2));
        PeriodicBoundaryXYZ<traitsT>::set_system(low, upp);
    }

    if(integral == "Underdamped Langevin")
    {
        return read_underdamped_langevin<traitsT>(sim, boundary, rng);
    }
    else if(integral == "NVE Newtonian")
    {
        return read_nve_newtonian<traitsT>(sim, boundary);
    }
    else
        throw std::invalid_argument("unknown time integral method: " + integral);
}


} // mjolnir
#endif /* MJOLNIR_TOML_READ_INTEGRATOR */
