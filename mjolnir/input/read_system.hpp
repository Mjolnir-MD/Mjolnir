#ifndef MJOLNIR_INPUT_READ_SYSTEM_HPP
#define MJOLNIR_INPUT_READ_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/math/vector_util.hpp>
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

    check_keys_available(system, {"boundary_shape"_s, "attributes"_s, "particles"_s});

    const auto& boundary  = toml::find<toml::value>(system, "boundary_shape");
    const auto& particles = toml::find<toml::array>(system, "particles");

    System<traitsT> sys(particles.size(), read_boundary<traitsT>(boundary));

    if(not particles.empty())
    {
        const auto& p = particles.front().as_table();
        // if the first particle has velocity, mjolnir assumes that velocity is
        // already initialized (velocity generation not required).
        sys.velocity_initialized() = (p.count("vel")      != 0) ||
                                     (p.count("velocity") != 0);
    }

    const auto find_either = [](const toml::value& v,
        const std::string& key1, const std::string& key2) -> const toml::value&
        {
            const auto& table = v.as_table();
            if(table.count(key1) != 0 && table.count(key2) != 0)
            {
                throw_exception<std::runtime_error>(toml::format_error("[error]"
                    " key duplicates.", table.at(key1), "here", table.at(key2),
                    "this conflicts with the above value definition"));
            }
            if(table.count(key1) != 0) {return table.at(key1);}
            if(table.count(key2) != 0) {return table.at(key2);}

            std::ostringstream oss;
            oss << "[error] both keys \"" << key1 << "\" and \"" << key2
                << "\" not found.";
            throw_exception<std::out_of_range>(
                    toml::format_error(oss.str(), v, "in this table"));
        };

    MJOLNIR_LOG_NOTICE(particles.size(), " particles are found.");
    for(std::size_t i=0; i<particles.size(); ++i)
    {
        using vec_type = std::array<real_type, 3>;
        const auto& p = particles.at(i);

        check_keys_available(p, {"mass"_s, "m"_s, "position"_s, "pos"_s,
                "velocity"_s, "vel"_s, "name"_s, "group"_s});

        sys.mass(i)     = toml::get<real_type>(find_either(p, "m",   "mass"));
        sys.position(i) = toml::get< vec_type>(find_either(p, "pos", "position"));
        sys.velocity(i) = math::make_coordinate<coordinate_type>(0, 0, 0);

        if(sys.velocity_initialized()) // it requires `velocity`.
        {
            sys.velocity(i) = toml::get<vec_type>(find_either(p, "vel", "velocity"));
        }
        else // velocity should not be there. check it.
        {
            if(p.as_table().count("vel") != 0 || p.as_table().count("velocity") != 0)
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
System<traitsT> read_system(const toml::value& root, std::size_t N)
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

        if(system_file.as_table().count("systems") != 1)
        {
            throw_exception<std::out_of_range>("[error] mjolnir::read_system: "
                "table [[systems]] not found in toml file\n --> ",
                input_path, file_name, "\n | the file should define [[systems]]"
                " table and define values in it.");
        }
        if(system_file.as_table().at("systems").is_array())
        {
            return read_system_from_table<traitsT>(
                toml::find<toml::array>(system_file, "systems").front());
        }
        return read_system_from_table<traitsT>(
                toml::find(system_file, "systems"));
    }
    return read_system_from_table<traitsT>(system_params.at(N));
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template System<SimulatorTraits<double, UnlimitedBoundary>       > read_system_from_table<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value& system);
extern template System<SimulatorTraits<float,  UnlimitedBoundary>       > read_system_from_table<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value& system);
extern template System<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_system_from_table<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& system);
extern template System<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_system_from_table<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& system);

extern template System<SimulatorTraits<double, UnlimitedBoundary>       > read_system<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<float,  UnlimitedBoundary>       > read_system<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_system<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_system<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& root, std::size_t N);
#endif

}//mjolnir
#endif //MJOLNIR_READ_SYSTEM
