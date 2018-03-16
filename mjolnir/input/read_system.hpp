#ifndef MJOLNIR_READ_SYSTEM
#define MJOLNIR_READ_SYSTEM
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <extlib/toml/toml.hpp>

namespace mjolnir
{

template<typename boundaryT> struct read_boundary_impl;

template<typename realT, typename coordT>
struct read_boundary_impl<UnlimitedBoundary<realT, coordT>>
{
    static UnlimitedBoundary<realT, coordT>
    invoke(const toml::Table& boundary)
    {
        return UnlimitedBoundary<realT, coordT>{};
    }
};

template<typename realT, typename coordT>
struct read_boundary_impl<CubicPeriodicBoundary<realT, coordT>>
{
    static CubicPeriodicBoundary<realT, coordT>
    invoke(const toml::Table& boundary)
    {
        const coordT upper(toml::get<std::array<realT, 3>>(
                    toml_value_at(boundary, "upper", "[boundary]")));
        const coordT lower(toml::get<std::array<realT, 3>>(
                    toml_value_at(boundary, "lower", "[boundary]")));
        return CubicPeriodicBoundary<realT, coordT>(lower, upper);
    }
};

template<typename traitsT>
typename traitsT::boundary_type
read_boundary(const toml::Table& boundary)
{
    return read_boundary_impl<typename traitsT::boundary_type>::invoke(boundary);
}

template<typename traitsT>
std::vector<typename System<traitsT>::particle_type>
read_particles(const toml::Table& system)
{
    typedef typename traitsT::real_type real_type;
    typedef typename traitsT::coordinate_type coordinate_type;

    std::vector<typename System<traitsT>::particle_type> ps;
    const auto& particles = toml_value_at(system, "particles", "[system]"
            ).cast<toml::value_t::Array>();
    ps.reserve(particles.size());

    for(const auto& p : particles)
    {
        const auto& params = p.cast<toml::value_t::Table>();
        const auto  mass = toml::get<real_type>(
                toml_value_at(params, "mass", "<anonymous> in particles"));
        const coordinate_type pos = toml::get<std::array<real_type, 3>>(
                toml_value_at(params, "position", "<anonymous> in particles"));
        const coordinate_type vel = toml::get<std::array<real_type, 3>>(
                toml_value_at(params, "velocity", "<anonymous> in particles"));
        const coordinate_type f(0.,0.,0.);
        ps.emplace_back(make_particle(mass, pos, vel, f));
    }
    return ps;
}

template<typename traitsT>
System<traitsT> read_system(const toml::Table& data, std::size_t N)
{
    typedef typename traitsT::real_type real_type;

    const auto& system_params = toml_value_at(data, "systems", "<root>"
            ).cast<toml::value_t::Array>();
    if(system_params.size() <= N)
    {
        throw std::out_of_range("no enough systems: " + std::to_string(N));
    }

    const auto& system = system_params.at(N).cast<toml::value_t::Table>();

    System<traitsT> sys(read_particles<traitsT>(system),
        read_boundary<traitsT>(toml_value_at(system, "boundary"
                ).cast<toml::value_t::Table>()));

    const auto& attributes = toml_value_at(system, "attributes", "[[systems]]"
            ).cast<toml::value_t::Table>();
    for(const auto& attr : attributes)
    {
        sys.attribute(attr.first) = toml::get<real_type>(attr.second);
    }
    return sys;
}

}//mjolnir
#endif //MJOLNIR_READ_SIMULATOR
