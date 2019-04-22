#ifndef MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/read_utility.hpp>
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
read_ignored_molecule(const toml::value& ignored_mol)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto name = toml::get<std::string>(ignored_mol);

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
            ignored_mol, "expected (Nothing|Self|Others)."));
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

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto radius_reader = [](const toml::value& v) -> real_type {
        return toml::find<real_type>(v, "radius");
    };
    auto parameters = read_array<real_type>(
            toml::find(global, "parameters"), radius_reader);

    return ExcludedVolumePotential<realT>(
        eps, std::move(parameters), ignore_particle_within,
        read_ignored_molecule(toml::find<toml::value>(ignore, "molecule")));
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

    using parameter_type = std::pair<real_type, real_type>;

    const auto lj_reader = [](const toml::value& v) -> parameter_type
    {
        const auto sigma   = toml::expect<real_type>(v, u8"σ").or_other(
                             toml::expect<real_type>(v, "sigma")).unwrap();
        const auto epsilon = toml::expect<real_type>(v, u8"ε").or_other(
                             toml::expect<real_type>(v, "epsilon")).unwrap();
        return std::make_pair(sigma, epsilon);
    };

    auto parameters = read_array<parameter_type>(
            toml::find(global, "parameters"), lj_reader);

    return LennardJonesPotential<realT>(
        std::move(parameters), ignore_particle_within,
        read_ignored_molecule(toml::find<toml::value>(ignore, "molecule")));
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

    const auto sigma   = toml::expect<real_type>(global, u8"σ").or_other(
                         toml::expect<real_type>(global, "sigma")).unwrap();
    const auto epsilon = toml::expect<real_type>(global, u8"ε").or_other(
                         toml::expect<real_type>(global, "epsilon")).unwrap();

    MJOLNIR_LOG_INFO("sigma   = ", sigma);
    MJOLNIR_LOG_INFO("epsilon = ", epsilon);

    return UniformLennardJonesPotential<realT>(
        sigma, epsilon, ignore_particle_within,
        read_ignored_molecule(toml::find<toml::value>(ignore, "molecule")));
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

    const auto charge_reader = [](const toml::value& v) -> real_type {
        return toml::find<real_type>(v, "charge");
    };
    auto parameters = read_array<real_type>(
            toml::find(global, "parameters"), charge_reader);

    return DebyeHuckelPotential<realT>(
        std::move(parameters), ignore_particle_within,
        read_ignored_molecule(toml::find<toml::value>(ignore, "molecule")));
}

} // mjolnir
#endif // MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
