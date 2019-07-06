#ifndef MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// ============================================================================
// global potential
// ============================================================================

inline IgnoreMolecule<typename Topology::molecule_id_type>
read_ignored_molecule(const toml::value& ignore)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto name = toml::find<std::string>(ignore, "molecule");

    if(name == "Nothing")
    {
        MJOLNIR_LOG_INFO("all the interactions"
                         "(both (inter|intra)-molecule) are included");
        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreNothing<typename Topology::molecule_id_type>>()
        };
    }
    else if(name == "Self" || name == "Intra")
    {
        MJOLNIR_LOG_INFO("intra-molecule interaction is ignored");
        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreSelf<typename Topology::molecule_id_type>>()
        };
    }
    else if(name == "Others" || name == "Inter")
    {
        MJOLNIR_LOG_INFO("inter-molecule interaction is ignored");
        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreOthers<typename Topology::molecule_id_type>>()
        };
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error(
            "[error] mjolnir::read_ignored_molecule: unknown setting",
            toml::find(ignore, "molecule"), "expected (Nothing|Self|Others)."));
    }
}

inline IgnoreGroup<typename Topology::group_id_type>
read_ignored_group(const toml::value& ignore)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using group_id_type = typename Topology::group_id_type;

    if(ignore.as_table().count("group") == 0)
    {
        MJOLNIR_LOG_NOTICE("No `ignore.group` is provided. "
                           "All the groups are taken into account.");
        return IgnoreGroup<group_id_type>{
            make_unique<IgnoreNoGroup<group_id_type>>()
        };
    }
    const auto name = toml::find<std::string>(ignore, "group");

    if(name == "Nothing")
    {
        MJOLNIR_LOG_INFO("all the interactions"
                         "(both (inter|intra)-group) are included");
        return IgnoreGroup<group_id_type>{
            make_unique<IgnoreNoGroup<group_id_type>>()
        };
    }
    else if(name == "Self" || name == "Intra")
    {
        MJOLNIR_LOG_INFO("intra-molecule interaction is ignored");
        return IgnoreGroup<group_id_type>{
            make_unique<IgnoreIntraGroup<group_id_type>>()
        };
    }
    else if(name == "Others" || name == "Inter")
    {
        MJOLNIR_LOG_INFO("inter-molecule interaction is ignored");
        return IgnoreGroup<group_id_type>{
            make_unique<IgnoreInterGroup<group_id_type>>()
        };
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error(
            "[error] mjolnir::read_ignored_group: unknown setting",
            toml::find(ignore, "gruop"),
            "expected (Nothing|(Self|Intra)|(Others|Inter))."));
    }
}

template<typename realT>
ExcludedVolumePotential<realT>
read_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using parameter_type = typename ExcludedVolumePotential<realT>::parameter_type;

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index");
        const auto radius = find_parameter<real_type  >(param, env, "radius");

        params.emplace_back(idx, radius);
        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", radius);
    }

    return ExcludedVolumePotential<realT>(
        eps, params, ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

template<typename realT>
LennardJonesPotential<realT>
read_lennard_jones_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using parameter_type = typename LennardJonesPotential<realT>::parameter_type;

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx     = find_parameter<std::size_t>(param, env, "index");
        const auto sigma   = find_parameter<real_type>(param, env, "sigma",   u8"σ");
        const auto epsilon = find_parameter<real_type>(param, env, "epsilon", u8"ε");

        params.emplace_back(idx, parameter_type{sigma, epsilon});
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma, ", epsilon = ", epsilon);
    }

    return LennardJonesPotential<realT>(
        std::move(params), ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

template<typename realT>
UniformLennardJonesPotential<realT>
read_uniform_lennard_jones_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto sigma   = toml::expect<real_type>(global, u8"σ").or_other(
                         toml::expect<real_type>(global, "sigma")).unwrap();
    const auto epsilon = toml::expect<real_type>(global, u8"ε").or_other(
                         toml::expect<real_type>(global, "epsilon")).unwrap();

    MJOLNIR_LOG_INFO("sigma   = ", sigma);
    MJOLNIR_LOG_INFO("epsilon = ", epsilon);

    using parameter_type = typename UniformLennardJonesPotential<realT>::parameter_type;
    std::vector<std::pair<std::size_t, parameter_type>> params;
    if(global.as_table().count("parameters") == 1)
    {
        for(const auto& param : toml::find<toml::array>(global, "parameters"))
        {
            const auto idx = find_parameter<std::size_t>(param, env, "index");
            params.emplace_back(idx, parameter_type{});
        }
    }
    return UniformLennardJonesPotential<realT>(
        sigma, epsilon, params, ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

template<typename realT>
DebyeHuckelPotential<realT>
read_debye_huckel_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using parameter_type = typename DebyeHuckelPotential<realT>::parameter_type;

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index");
        const auto charge = find_parameter<real_type  >(param, env, "charge");

        params.emplace_back(idx, parameter_type{charge});
        MJOLNIR_LOG_INFO("idx = ", idx, ", charge = ", charge);
    }
    return DebyeHuckelPotential<realT>(
        std::move(params), ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

} // mjolnir
#endif // MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
