#ifndef MJOLNIR_READ_SYSTEM
#define MJOLNIR_READ_SYSTEM
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/get_toml_value.hpp>

namespace mjolnir
{
namespace detail
{

template<typename boundaryT> struct read_boundary_impl;

template<typename realT, typename coordT>
struct read_boundary_impl<UnlimitedBoundary<realT, coordT>>
{
    static UnlimitedBoundary<realT, coordT>
    invoke(const toml::Table& boundary)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(read_boundary_impl::invoke(), 0);
        MJOLNIR_LOG_INFO("no boundary is set");
        return UnlimitedBoundary<realT, coordT>{};
    }
};

template<typename realT, typename coordT>
struct read_boundary_impl<CubicPeriodicBoundary<realT, coordT>>
{
    static CubicPeriodicBoundary<realT, coordT>
    invoke(const toml::Table& boundary)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(read_boundary_impl::invoke(), 0);
        MJOLNIR_LOG_INFO("shape of periodic boundary is cubic");

        const coordT upper(toml::get<std::array<realT, 3>>(
                    toml_value_at(boundary, "upper", "[boundary]")));
        const coordT lower(toml::get<std::array<realT, 3>>(
                    toml_value_at(boundary, "lower", "[boundary]")));
        assert(upper[0] > lower[0]);
        assert(upper[1] > lower[1]);
        assert(upper[2] > lower[2]);
        MJOLNIR_LOG_INFO("upper limit of the boundary = ", upper);
        MJOLNIR_LOG_INFO("lower limit of the boundary = ", lower);
        return CubicPeriodicBoundary<realT, coordT>(lower, upper);
    }
};

} // detail

template<typename traitsT>
typename traitsT::boundary_type
read_boundary(const toml::Table& boundary)
{
    using boundary_t = typename traitsT::boundary_type;
    return detail::read_boundary_impl<boundary_t>::invoke(boundary);
}

template<typename traitsT>
System<traitsT> read_system(const toml::Table& data, std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_system(), 0);

    using real_type  = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const auto& system_params = toml_value_at(data, "systems", "<root>"
            ).cast<toml::value_t::Array>();
    if(system_params.size() <= N)
    {
        throw_exception<std::out_of_range>("no enough system definitions: ", N);
    }

    MJOLNIR_LOG_INFO(system_params.size(), " systems are provided");
    MJOLNIR_LOG_INFO("using ", N, "-th system");

    const auto& system   = system_params.at(N).cast<toml::value_t::Table>();
    const auto& boundary = toml_value_at(system, "boundary","[systems.boundary]"
            ).cast<toml::value_t::Table>();
    const auto& particles = toml_value_at(system, "particles", "[system]"
            ).cast<toml::value_t::Array>();

    System<traitsT> sys(particles.size(), read_boundary<traitsT>(boundary));
    MJOLNIR_LOG_INFO(particles.size(), " particles are provided");
    for(std::size_t i=0; i<particles.size(); ++i)
    {
        using vec_type = std::array<real_type, 3>;
        const auto& params = particles[i].cast<toml::value_t::Table>();

        sys[i].mass     = toml::get<real_type>(toml_value_at(params, "mass", "element of [[system.particles]]"));
        sys[i].position = toml::get< vec_type>(toml_value_at(params, "position", "element of [[system.particles]]"));
        sys[i].velocity = toml::get< vec_type>(toml_value_at(params, "velocity", "element of [[system.particles]]"));
        sys[i].force    = coordinate_type(0, 0, 0);

        MJOLNIR_LOG_DEBUG("mass     = ", sys[i].mass    );
        MJOLNIR_LOG_DEBUG("position = ", sys[i].position);
        MJOLNIR_LOG_DEBUG("velocity = ", sys[i].velocity);
        MJOLNIR_LOG_DEBUG("force    = ", sys[i].force   );
    }

    const auto& attributes = toml_value_at(system, "attributes", "[[systems]]"
            ).cast<toml::value_t::Table>();
    for(const auto& attr : attributes)
    {
        const real_type attribute = toml::get<real_type>(attr.second);
        MJOLNIR_LOG_INFO("attribute.", attr.first, " = ", attribute);
        sys.attribute(attr.first) = attribute;
    }
    return sys;
}

}//mjolnir
#endif //MJOLNIR_READ_SIMULATOR
