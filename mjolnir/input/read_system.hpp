#ifndef MJOLNIR_INPUT_READ_SYSTEM_HPP
#define MJOLNIR_INPUT_READ_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/format_nth.hpp>
#include <mjolnir/math/vector_util.hpp>
#include <mjolnir/input/read_table_from_file.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/input/utility.hpp>

namespace mjolnir
{
namespace detail
{

// reads boundary settings depending on the type information.
// the specific cases are implemented below.
template<typename boundaryT> struct read_boundary_impl;

template<typename realT, typename coordT>
struct read_boundary_impl<UnlimitedBoundary<realT, coordT>>
{
    static UnlimitedBoundary<realT, coordT>
    invoke(const toml::value& v)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("no boundary is set. unlimited");
        check_keys_available(v, {}); // no available keys
        return UnlimitedBoundary<realT, coordT>{};
    }
};

template<typename realT, typename coordT>
struct read_boundary_impl<CuboidalPeriodicBoundary<realT, coordT>>
{
    static CuboidalPeriodicBoundary<realT, coordT>
    invoke(const toml::value& boundary)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("shape of periodic boundary is cuboid");
        check_keys_available(boundary, {"upper"_s, "lower"_s});

        const auto upper = toml::find<coordT>(boundary, "upper");
        const auto lower = toml::find<coordT>(boundary, "lower");

        if(math::X(upper) <= math::X(lower) ||
           math::Y(upper) <= math::Y(lower) || math::Z(upper) <= math::Z(lower))
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_boundary: upper should be larger than lower",
                toml::find(boundary, "upper"), "upper boundary here",
                toml::find(boundary, "lower"), "lower boundary here"));
        }
        MJOLNIR_LOG_INFO("upper limit of the boundary = ", upper);
        MJOLNIR_LOG_INFO("lower limit of the boundary = ", lower);

        return CuboidalPeriodicBoundary<realT, coordT>(lower, upper);
    }
};

} // detail

template<typename traitsT>
typename traitsT::boundary_type
read_boundary(const toml::value& boundary)
{
    using boundary_t = typename traitsT::boundary_type;
    return detail::read_boundary_impl<boundary_t>::invoke(boundary);
}

// It reads particles and other system-specific attributes (e.g. temperature)
template<typename traitsT>
System<traitsT>
read_system(const toml::value& root, const std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const auto& systems = toml::find(root, "systems");

    MJOLNIR_LOG_NOTICE("reading ", format_nth(N), " system ...");
    const auto system = read_table_from_file(
            systems.at(N), "systems", read_input_path(root));

    check_keys_available(system, {"boundary_shape"_s, "attributes"_s, "particles"_s});

    const auto& boundary  = toml::find<toml::value>(system, "boundary_shape");
    const auto& particles = toml::find<toml::array>(system, "particles");

    System<traitsT> sys(particles.size(), read_boundary<traitsT>(boundary));

    for(const auto& attr : toml::find<toml::table>(system, "attributes"))
    {
        const real_type attribute = toml::get<real_type>(attr.second);
        sys.attribute(attr.first) = attribute;
        MJOLNIR_LOG_INFO("attribute.", attr.first, " = ", attribute);
    }

    if(particles.empty())
    {
        return sys;
    }

    {
        // if the first particle has velocity, mjolnir assumes that velocity is
        // already initialized (velocity generation not required).
        const auto& p = particles.front().as_table();
        sys.velocity_initialized() = (p.count("vel")      != 0) ||
                                     (p.count("velocity") != 0);
    }

    const auto find_either = [](const toml::value& v, const std::string& key1,
                                const std::string& key2) -> const toml::value&
        {
            if(v.contains(key1) && v.contains(key2) != 0)
            {
                throw_exception<std::runtime_error>(toml::format_error("[error]"
                    " key duplicates.", v.at(key1), "here", v.at(key2),
                    "this conflicts with the above value definition"));
            }
            if(v.contains(key1)) {return v.at(key1);}
            if(v.contains(key2)) {return v.at(key2);}

            throw_exception<std::out_of_range>(toml::format_error(
                "both keys, \""_s + key1 + "\" and \""_s + key2 +
                "\", are not found."_s, v, "in this table"_s));
        };

    MJOLNIR_LOG_NOTICE(particles.size(), " particles are found.");
    for(std::size_t i=0; i<particles.size(); ++i)
    {
        using coord_type = coordinate_type;
        const auto& p = particles.at(i);

        check_keys_available(p, {"mass"_s, "m"_s, "position"_s, "pos"_s,
                "velocity"_s, "vel"_s, "name"_s, "group"_s});

        sys.mass(i)     = toml::get<real_type >(find_either(p, "m", "mass"));
        sys.position(i) = toml::get<coord_type>(find_either(p, "pos", "position"));
        sys.velocity(i) = math::make_coordinate<coordinate_type>(0, 0, 0);

        if(sys.velocity_initialized()) // it requires `velocity`.
        {
            sys.velocity(i) = toml::get<coord_type>(find_either(p, "vel", "velocity"));
        }
        else // if velocity_initialized is false, velocity should not be there.
        {
            if(p.contains("vel") || p.contains("velocity"))
            {
                throw_exception<std::runtime_error>(toml::format_error("[error]"
                    "read_system(): missing (or extraneous) velocity input",
                    particles.front(), "this does not have `velocity` field",
                    p, "but this has `velocity` field", {
                    "To specify the initial velocity, it is necessary to give "
                    "the value at all particles.",
                    "To generate the initial velocity with a random number, "
                    "all the `velocity` fields must be empty."
                    }));
            }
        }

        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        sys.name(i)  = toml::find_or<std::string>(p, "name",  "X"   );
        sys.group(i) = toml::find_or<std::string>(p, "group", "NONE");
        sys.rmass(i) = 1.0 / sys.mass(i);

        MJOLNIR_LOG_INFO("mass = ",        sys.mass(i),
                          ", position = ", sys.position(i),
                          ", velocity = ", sys.velocity(i),
                          ", force = ",    sys.force(i),
                          ", name = ",     sys.name(i),
                          ", group = ",    sys.group(i));
    }
    return sys;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template System<SimulatorTraits<double, UnlimitedBoundary>       > read_system<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<float,  UnlimitedBoundary>       > read_system<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_system<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_system<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& root, std::size_t N);
#endif

}//mjolnir
#endif //MJOLNIR_READ_SYSTEM
