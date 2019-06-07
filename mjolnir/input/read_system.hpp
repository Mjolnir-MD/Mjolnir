#ifndef MJOLNIR_INPUT_READ_SYSTEM_HPP
#define MJOLNIR_INPUT_READ_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/math/vector_util.hpp>
#include <mjolnir/input/read_files_table.hpp>

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
    invoke(const toml::value&)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("no boundary is set. unlimited");
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

        const auto upper = toml::find<std::array<realT, 3>>(boundary, "upper");
        const auto lower = toml::find<std::array<realT, 3>>(boundary, "lower");
        if(upper[0] <= lower[0] || upper[1] <= lower[1] || upper[2] <= lower[2])
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
System<traitsT> read_system_from_table(const toml::value& system)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    MJOLNIR_LOG_NOTICE("reading system ...");

    const auto& boundary  = toml::find<toml::value>(system, "boundary_shape");
    const auto& particles = toml::find<toml::array>(system, "particles");

    System<traitsT> sys(particles.size(), read_boundary<traitsT>(boundary));

    MJOLNIR_LOG_NOTICE(particles.size(), " particles are found.");
    for(std::size_t i=0; i<particles.size(); ++i)
    {
        using vec_type = std::array<real_type, 3>;
        const auto& p = particles.at(i);

        sys[i].mass     = toml::expect<real_type>(p, "m").or_other(
                          toml::expect<real_type>(p, "mass")).unwrap();
        sys[i].rmass    = 1.0 / sys[i].mass;
        sys[i].position = toml::expect<vec_type>(p, "pos").or_other(
                          toml::expect<vec_type>(p, "position")).unwrap();
        sys[i].velocity = toml::expect<vec_type>(p, "vel").or_other(
                          toml::expect<vec_type>(p, "velocity")).unwrap();
        sys[i].force    = math::make_coordinate<coordinate_type>(0, 0, 0);
        sys[i].name     = toml::expect<std::string>(p, "name").unwrap_or("X");
        sys[i].group    = toml::expect<std::string>(p, "group").unwrap_or("NONE");

        MJOLNIR_LOG_INFO("mass = ",        sys[i].mass,
                          ", position = ", sys[i].position,
                          ", velocity = ", sys[i].velocity,
                          ", force = ",    sys[i].force,
                          ", name = ",     sys[i].name,
                          ", group = ",    sys[i].group);
    }
    MJOLNIR_LOG_NOTICE("done.");

    for(const auto& attr : toml::find<toml::table>(system, "attributes"))
    {
        const real_type attribute = toml::get<real_type>(attr.second);
        sys.attribute(attr.first) = attribute;
        MJOLNIR_LOG_INFO("attribute.", attr.first, " = ", attribute);
    }
    return sys;
}

// reads N-th system. In some (abnormal) simulation, e.g. REMD or string method,
// we may have N different systems. This function reads N-th system. but in most
// of the cases, N is always equal to 0.
template<typename traitsT>
System<traitsT> read_system(const toml::table& root, std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    if(N != 0) {MJOLNIR_LOG_NOTICE("reading ", N, "-th [[system]].");}
    else       {MJOLNIR_LOG_NOTICE("reading [[system]].");}

    const auto& system_params = toml::find<toml::array>(root, "systems");
    if(system_params.size() <= N)
    {
        throw_exception<std::out_of_range>("[error] mjolnir::read_system: "
            "no enough system definitions: ", N, " is required, but only ",
            system_params.size(), " defined.");
    }
    MJOLNIR_LOG_INFO(system_params.size(), " systems are provided");
    MJOLNIR_LOG_INFO("using ", N, "-th system");

    const auto& system = toml::get<toml::table>(system_params.at(N));
    if(system.count("file_name") == 1)
    {
        const auto input_path = read_input_path(root);
        const auto file_name = toml::get<std::string>(system.at("file_name"));
        MJOLNIR_LOG_NOTICE("system is defined in ", input_path + file_name);
        if(system.size() != 1)
        {
            MJOLNIR_LOG_WARN("[[systems]] has \"file_name\" and other values.");
            MJOLNIR_LOG_WARN("When \"file_name\" is provided, other values are "
                             "ignored because those are read from the specified"
                             " file (", input_path, file_name, ").");
        }

        MJOLNIR_LOG_NOTICE("reading ", input_path, file_name, " ...");
        const auto system_file = toml::parse(input_path + file_name);
        MJOLNIR_LOG_NOTICE(" done.");

        if(system_file.count("systems") != 1)
        {
            throw_exception<std::out_of_range>("[error] mjolnir::read_system: "
                "table [[systems]] not found in toml file\n --> ",
                input_path, file_name, "\n | the file should define [[systems]]"
                " table and define values in it.");
        }
        if(system_file.at("systems").is(toml::value_t::Array))
        {
            return read_system_from_table<traitsT>(
                toml::find<toml::array>(system_file, "systems").front());
        }
        return read_system_from_table<traitsT>(
                toml::find(system_file, "systems"));
    }
    return read_system_from_table<traitsT>(system_params.at(N));
}

}//mjolnir
#endif //MJOLNIR_READ_SYSTEM
