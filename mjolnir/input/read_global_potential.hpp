#ifndef MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
#define MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// ============================================================================
// global potential
// ============================================================================

inline IgnoreMolecule<typename Topology::molecule_id_type>
read_ignored_molecule(const std::string& ignored_mol)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_ignored_molecule(), 0);

    if(ignored_mol == "Nothing")
    {
        MJOLNIR_LOG_INFO("all the interactions(both (inter|intra)-molecule) are included");
        return make_unique<IgnoreNothing<typename Topology::molecule_id_type>>();
    }
    else if(ignored_mol == "Self" || ignored_mol == "Intra")
    {
        MJOLNIR_LOG_INFO("intra-molecule interaction is ignored");
        return make_unique<IgnoreSelf<typename Topology::molecule_id_type>>();
    }
    else if(ignored_mol == "Others" || ignored_mol == "Inter")
    {
        MJOLNIR_LOG_INFO("inter-molecule interaction is ignored");
        return make_unique<IgnoreOthers<typename Topology::molecule_id_type>>();
    }
    else
    {
        throw_exception<std::runtime_error>("invalid `ignored_molecule`: ",
            ignored_mol, ". allowed: Nothing, Self, or Others.");
    }
}

template<typename realT>
ExcludedVolumePotential<realT>
read_excluded_volume_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_excluded_volume_potential(), 0);
    using real_type = realT;

    const auto location = "[forcefield.global] for ExcludedVolume potential";

    const auto& ignore =
        get_toml_value<toml::Table>(global, "ignore", location);

    auto ignored_mol = read_ignored_molecule(
        get_toml_value<std::string>(ignore, "molecule", location));

    std::map<std::string, std::size_t> connections;
    for(const auto connection :
            get_toml_value<toml::Table>(ignore, "particles_within", location))
    {
        connections[connection.first] =
            toml::get<std::size_t>(connection.second);
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
                         " within ", connections.at(connection.first), " will ",
                         "be ignored");
    }

    const real_type eps = get_toml_value<real_type>(
        global, "epsilon", location);
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& ps = get_toml_value<toml::Array>(global, "parameters", location);
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");
    const auto parameters =
        "element of [[forcefield.global.parameters]] for ExcludedVolume";

    std::vector<real_type> params;
    params.reserve(ps.size());

    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto  idx = get_toml_value<std::size_t>(tab, "index", parameters);
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }

        const auto radius = get_toml_value<real_type>(tab, "radius", parameters);
        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", radius);
        params.at(idx) = radius;
    }

    return ExcludedVolumePotential<realT>(
        eps, std::move(params), connections, std::move(ignored_mol));
}

template<typename realT>
LennardJonesPotential<realT>
read_lennard_jones_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_lennard_jones_potential(), 0);
    using real_type = realT;

    const auto location   = "[forcefield.global] for LennardJones potential";

    const auto& ignore =
        get_toml_value<toml::Table>(global, "ignore", location);

    auto ignored_mol = read_ignored_molecule(
        get_toml_value<std::string>(ignore, "molecule", location));

    std::map<std::string, std::size_t> connections;
    for(const auto connection :
        get_toml_value<toml::Table>(ignore, "particles_within", location))
    {
        connections[connection.first] =
            toml::get<std::size_t>(connection.second);
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
                         " within ", connections.at(connection.first), " will ",
                         "be ignored");
    }

    const auto& ps = get_toml_value<toml::Array>(global, "parameters", location);
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");
    const auto parameters =
        "element of [[forcefield.global.parameters]] for Lennard-Jones";

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index", parameters);
        if(params.size() <= idx)
        {
            const std::pair<real_type, real_type> dummy{0., 0.};
            params.resize(idx+1, dummy);
        }
        const auto sigma   = get_toml_value<real_type>(tab, {"sigma"_s,   u8"σ"_s}, parameters);
        const auto epsilon = get_toml_value<real_type>(tab, {"epsilon"_s, u8"ε"_s}, parameters);
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma, ", epsilon = ", epsilon);

        params.at(idx) = std::make_pair(sigma, epsilon);
    }

    return LennardJonesPotential<realT>(
            std::move(params), connections, std::move(ignored_mol));
}

template<typename realT>
UniformLennardJonesPotential<realT>
read_uniform_lennard_jones_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_uniform_lennard_jones_potential(), 0);
    using real_type = realT;
    const auto location = "[forcefield.global] for UniformLennardJones";

    const auto& ignore =
        get_toml_value<toml::Table>(global, "ignore", location);

    auto ignored_mol = read_ignored_molecule(
        get_toml_value<std::string>(ignore, "molecule", location));

    std::map<std::string, std::size_t> connections;
    for(const auto connection :
            get_toml_value<toml::Table>(ignore, "particles_within", location))
    {
        connections[connection.first] =
            toml::get<std::size_t>(connection.second);
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
                         " within ", connections.at(connection.first), " will ",
                         "be ignored");
    }

    const auto sigma   = get_toml_value<real_type>(global, {"sigma"_s,   u8"σ"_s}, location);
    const auto epsilon = get_toml_value<real_type>(global, {"epsilon"_s, u8"ε"_s}, location);

    MJOLNIR_LOG_INFO("sigma   = ", sigma);
    MJOLNIR_LOG_INFO("epsilon = ", epsilon);

    return UniformLennardJonesPotential<realT>(
            sigma, epsilon, connections, std::move(ignored_mol));
}

template<typename realT>
DebyeHuckelPotential<realT>
read_debye_huckel_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_debye_huckel_potential(), 0);
    using real_type = realT;
    const auto location = "[forcefield.global] for DebyeHuckel";

    const auto& ignore =
        get_toml_value<toml::Table>(global, "ignore", location);

    auto ignored_mol = read_ignored_molecule(
        get_toml_value<std::string>(ignore, "molecule", location));

    std::map<std::string, std::size_t> connections;
    for(const auto connection :
            get_toml_value<toml::Table>(ignore, "particles_within", location))
    {
        connections[connection.first] =
            toml::get<std::size_t>(connection.second);
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
                         " within ", connections.at(connection.first), " will ",
                         "be ignored");
    }

    const auto& ps = get_toml_value<toml::Array>(global, "parameters", location);
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    const auto parameters =
        "element of [[forcefield.global.parameters]] for Debye-Huckel";

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index", parameters);
        const auto charge = get_toml_value<real_type>(tab, "charge", parameters);
        MJOLNIR_LOG_INFO("idx    = ", idx, ", charge = ", charge);
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        params.at(idx) = charge;
    }

    return DebyeHuckelPotential<realT>(
            std::move(params), connections, std::move(ignored_mol));
}

} // mjolnir
#endif // MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
