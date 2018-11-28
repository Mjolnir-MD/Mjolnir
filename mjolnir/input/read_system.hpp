#ifndef MJOLNIR_READ_SYSTEM
#define MJOLNIR_READ_SYSTEM
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_files_table.hpp>

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
        MJOLNIR_LOG_INFO("no boundary is set. unlimited");
        return UnlimitedBoundary<realT, coordT>{};
    }
};

template<typename realT, typename coordT>
struct read_boundary_impl<CuboidalPeriodicBoundary<realT, coordT>>
{
    static CuboidalPeriodicBoundary<realT, coordT>
    invoke(const toml::Table& boundary)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(read_boundary_impl::invoke(), 0);
        MJOLNIR_LOG_INFO("shape of periodic boundary is cuboid");

        const coordT upper(get_toml_value<std::array<realT, 3>>(
                    boundary, "upper", "[boundary_shape]"));
        const coordT lower(get_toml_value<std::array<realT, 3>>(
                    boundary, "lower", "[boundary_shape]"));

        assert(upper[0] > lower[0]);
        assert(upper[1] > lower[1]);
        assert(upper[2] > lower[2]);
        MJOLNIR_LOG_INFO("upper limit of the boundary = ", upper);
        MJOLNIR_LOG_INFO("lower limit of the boundary = ", lower);
        return CuboidalPeriodicBoundary<realT, coordT>(lower, upper);
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
System<traitsT> read_system_from_table(const toml::Table& system)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_system_from_table(), 0);
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const auto& boundary  =
        get_toml_value<toml::Table>(system, "boundary_shape", "[system]");
    const auto& particles =
        get_toml_value<toml::Array>(system, "particles", "[system]");

    System<traitsT> sys(particles.size(), read_boundary<traitsT>(boundary));
    MJOLNIR_LOG_NOTICE(particles.size(), " particles are found. reading...");
    for(std::size_t i=0; i<particles.size(); ++i)
    {
        using vec_type = std::array<real_type, 3>;
        const auto& params = particles[i].cast<toml::value_t::Table>();

        sys[i].mass     = get_toml_value<real_type>(
                params, {"m","mass"},       "[[system.particles]]");
        sys[i].rmass    = 1.0 / sys[i].mass;
        sys[i].position = get_toml_value<vec_type>(
                params, {"pos","position"}, "[[system.particles]]");
        sys[i].velocity = get_toml_value<vec_type>(
                params, {"vel","velocity"}, "[[system.particles]]");
        sys[i].force    = coordinate_type(0, 0, 0);
        sys[i].name     = "X";
        sys[i].group    = "NONE";

        if(params.count("name") == 1)
        {
            sys[i].name = toml::get<std::string>(params.at("name"));
        }
        if(params.count("group") == 1)
        {
            sys[i].group = toml::get<std::string>(params.at("group"));
        }

        MJOLNIR_LOG_INFO("mass = ",        sys[i].mass,
                          ", position = ", sys[i].position,
                          ", velocity = ", sys[i].velocity,
                          ", force = ",    sys[i].force,
                          ", name = ",     sys[i].name,
                          ", group = ",    sys[i].group);
    }
    MJOLNIR_LOG_NOTICE("done.");

    const auto& attributes =
        get_toml_value<toml::Table>(system, "attributes", "[[systems]]");
    for(const auto& attr : attributes)
    {
        if(attr.second.type() != toml::value_t::Float)
        {
            throw_exception<std::runtime_error>("mjolnir::read_system: "
                "attribute `", attr.first, "` has type `", attr.second.type(),
                "`, expected type is `toml::Float (a.k.a. double)`.");
        }
        const real_type attribute = toml::get<real_type>(attr.second);
        MJOLNIR_LOG_INFO("attribute.", attr.first, " = ", attribute);
        sys.attribute(attr.first) = attribute;
    }
    return sys;
}

template<typename traitsT>
System<traitsT> read_system(const toml::Table& data, std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_system(), 0);
    MJOLNIR_LOG_NOTICE("reading ", N, "-th [[system]].");

    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const auto& system_params =
        get_toml_value<toml::Array>(data, "systems", "<root>");
    if(system_params.size() <= N)
    {
        throw_exception<std::out_of_range>("no enough system definitions: ", N);
    }

    const auto& files = get_toml_value<toml::Table>(data, "files", "<root>");
    const auto input_path = read_input_path(data);

    MJOLNIR_LOG_INFO(system_params.size(), " systems are provided");
    MJOLNIR_LOG_INFO("using ", N, "-th system");

    const auto& system = system_params.at(N).cast<toml::value_t::Table>();
    if(system.count("file_name") == 1)
    {
        const auto file_name =
            get_toml_value<std::string>(system, "file_name", "[[systems]]");
        MJOLNIR_LOG_INFO("file_name = ", file_name);

        if(system.size() != 1)
        {
            MJOLNIR_LOG_WARN("[[systems]] has `file_name` key and other keys.");
            MJOLNIR_LOG_WARN("When `file_name` is provided, other values are "
                             "ignored because those are read from the specified"
                             " file (", file_name, ").");
        }

        MJOLNIR_LOG_NOTICE("system is defined in ", input_path + file_name);
        const auto system_file = toml::parse(input_path + file_name);
        if(system_file.count("systems") == 1)
        {
            MJOLNIR_LOG_WARN("in `system` file, root object is treated as "
                             "one of the [systems] tables.");
            MJOLNIR_LOG_WARN("but in ", file_name, ", [systems] table found."
                             "trying to read it as a system setup.");

            if(system_file.at("systems").type() != toml::value_t::Table)
            {
                MJOLNIR_LOG_ERROR("type of `systems` is different from "
                                  "toml::Table in file (", file_name, ").");
                MJOLNIR_LOG_ERROR("note: [[...]] means Array-of-Tables. "
                                  "please take care.");
                std::exit(1);
            }
            return read_system_from_table<traitsT>(get_toml_value<toml::Table>(
                    system_file, "systems", file_name));
        }
        return read_system_from_table<traitsT>(system_file);
    }
    return read_system_from_table<traitsT>(system);
}

}//mjolnir
#endif //MJOLNIR_READ_SYSTEM
