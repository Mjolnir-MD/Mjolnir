#ifndef MJOLNIR_TOML_READ_PARTICLES
#define MJOLNIR_TOML_READ_PARTICLES
#include <mjolnir/core/ParticleContainer.hpp>
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
ParticleContainer<traitsT>
read_particles(const toml::Table& tab)
{
    typedef typename traitsT::coordinate_type coordT;

    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_particles CALLED");

    const std::size_t num_particles = static_cast<std::size_t>(
        toml::get<toml::Integer>(tab.at("number_of_particle")));
    MJOLNIR_LOG_INFO("num_particles", num_particles);

    const std::vector<toml::Table> particles =
        toml::get<toml::Array<toml::Table>>(tab.at("particles"));
    MJOLNIR_LOG_INFO("particles.size()", particles.size());
    if(num_particles != particles.size())
        throw std::runtime_error(
                "number_of_particles is not equal to particles.size()");

    ParticleContainer<traitsT> pcon(particles.size());

    for(auto iter = make_zip(particles.cbegin(), pcon.begin());
            iter != make_zip(particles.cend(), pcon.end()); ++iter)
    {
        const auto p = get<0>(iter);
        const auto mass = toml::get<toml::Float>(p->at("mass"));
        const auto pos  = toml::get<toml::Array<toml::Float>>(p->at("position"));
        const auto vel  = toml::get<toml::Array<toml::Float>>(p->at("velocity"));

        get<1>(iter)->mass     = mass;
        get<1>(iter)->position = coordT(pos.at(0), pos.at(1), pos.at(2));
        get<1>(iter)->velocity = coordT(vel.at(0), vel.at(1), vel.at(2));
        get<1>(iter)->force    = coordT(0., 0., 0.);
    }
    MJOLNIR_LOG_DEBUG("read_particles RETURNED");
    return pcon;
}




} // mjolnir
#endif /* MJOLNIR_TOML_READ_PARTICLES */
