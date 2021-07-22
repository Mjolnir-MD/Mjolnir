#ifndef MJOLNIR_INPUT_READ_SYSTEM_HPP
#define MJOLNIR_INPUT_READ_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/MsgPackLoader.hpp>
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
    invoke(const toml::value& system)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("no boundary is set. unlimited");

        if(system.contains("boundary_shape"))
        {
            check_keys_available(system.at("boundary_shape"), {/*no keys*/});
        }
        return UnlimitedBoundary<realT, coordT>{};
    }
};

template<typename realT, typename coordT>
struct read_boundary_impl<CuboidalPeriodicBoundary<realT, coordT>>
{
    static CuboidalPeriodicBoundary<realT, coordT>
    invoke(const toml::value& system)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("shape of periodic boundary is cuboid");

        const auto& boundary = system.at("boundary_shape");
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
read_boundary(const toml::value& system)
{
    using boundary_t = typename traitsT::boundary_type;
    return detail::read_boundary_impl<boundary_t>::invoke(system);
}

template<typename traitsT>
System<traitsT> load_system_from_msgpack(const std::string& msg_file)
{
    MsgPackLoader<traitsT> loader;
    return loader.load_system(msg_file);
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

    // ------------------------------------------------------------------------
    // check [[systems]] has `file_name = "checkpoint.msg"`.
    // If `file_name` has `.msg` file, load system status from the file.
    if(systems.at(N).contains("file_name"))
    {
        if(1 < systems.at(N).size())
        {
            MJOLNIR_LOG_WARN("When `file_name` is specified, "
                             "all other keys are ignored.");
            check_keys_available(systems.at(N), {"file_name"_s});
        }

        const auto& fname = toml::find<std::string>(systems, N, "file_name");
        if(file_extension_is(fname, ".msg"))
        {
            const auto& input_path = get_input_path();

            MJOLNIR_LOG_NOTICE("msgpack file specified. load system status from ",
                               input_path, fname);

            return load_system_from_msgpack<traitsT>(input_path + fname);
        }
        else
        {
            MJOLNIR_LOG_ERROR("unknown file format: ", fname, ".");
            MJOLNIR_LOG_ERROR("supported format is: .msg");
            throw_exception<std::runtime_error>("[error] mjolnir::read_system: "
                    "mjolnir supports .msg file for restarting");
        }
    }

    // ------------------------------------------------------------------------
    // `[[systems]]` does not have "checkpoint.msg". It should have all the
    // settings. read it.

    const auto system = read_table_from_file(systems.at(N), "systems");

    check_keys_available(system, {"boundary_shape"_s, "attributes"_s, "dynamic_variables"_s, "particles"_s});

    const auto& particles = toml::find<toml::array>(system, "particles");

    System<traitsT> sys(particles.size(), read_boundary<traitsT>(system));

    // -----------------------------------------------------------------------
    // read attributes and dynamic variables

    for(const auto& attr : toml::find<toml::table>(system, "attributes"))
    {
        const real_type attribute = toml::get<real_type>(attr.second);
        sys.attribute(attr.first) = attribute;
        MJOLNIR_LOG_INFO("attribute.", attr.first, " = ", attribute);
    }

    if(system.contains("dynamic_variables"))
    {
        for(const auto& dynvar : toml::find<toml::table>(system, "dynamic_variables"))
        {
            MJOLNIR_LOG_INFO("a dynamic varaible ", dynvar.first, " found");
            const auto& var = dynvar.second;

            const auto x = toml::find   <real_type>(var, "x");
            const auto m = toml::find   <real_type>(var, "m");
            const auto v = toml::find_or<real_type>(var, "v", real_type(0.0));
            const auto f = toml::find_or<real_type>(var, "f", real_type(0.0));
            const auto g = toml::find_or<real_type>(var, "gamma", real_type(0.0));

            if(var.contains("boundary"))
            {
                const auto& boundary = toml::find<std::string>(var, "boundary");
                const auto lower = toml::find<real_type>(var, "lower");
                const auto upper = toml::find<real_type>(var, "upper");

                if(boundary == "Periodic")
                {
                    MJOLNIR_LOG_INFO("boundary is periodic");
                    sys.variable(dynvar.first) = DynamicVariable<real_type>(
                        PeriodicDynamicVariable<real_type>(x, v, f, m, g, lower, upper));
                }
                else if(boundary == "Repulsive")
                {
                    MJOLNIR_LOG_INFO("boundary is repulsive");
                    sys.variable(dynvar.first) = DynamicVariable<real_type>(
                        RepulsiveDynamicVariable<real_type>(x, v, f, m, g, lower, upper));
                }
                else
                {
                    MJOLNIR_LOG_ERROR("unknown boundary setting: ", boundary);
                    throw_exception<std::runtime_error>("[error] mjolnir::read_system: "
                        "available boundaries of dynvar are: \"Periodic\" or \"Repulsive\"");
                }
            }
            else
            {
                MJOLNIR_LOG_INFO("no boundary");
                sys.variable(dynvar.first) = DynamicVariable<real_type>(
                            DefaultDynamicVariable<real_type>(x, v, f, m, g));
            }
            MJOLNIR_LOG_INFO("dynvar.", dynvar.first, ".m = ", m, " .x = ", x, " .v = ", v, " .f = ", f);
        }
    }

    // if there is no particle, return.
    if(particles.empty())
    {
        return sys;
    }

    // -----------------------------------------------------------------------
    // prepare for reading particles ...

    {
        // if the first particle has velocity, mjolnir assumes that velocity is
        // already initialized (velocity generation not required).
        const auto& p = particles.front();
        sys.velocity_initialized() = (p.contains("vel") || p.contains("velocity"));
    }

    // A functor to find a value that corresponds to either of the key.
    // If both key exists, throw an error.
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

    // -----------------------------------------------------------------------
    // read all the particles ...

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
