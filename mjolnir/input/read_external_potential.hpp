#ifndef MJOLNIR_READ_EXTERNAL_POTENTIAL_HPP
#define MJOLNIR_READ_EXTERNAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/external/ImplicitMembranePotential.hpp>
#include <mjolnir/potential/external/LennardJonesWallPotential.hpp>
#include <mjolnir/potential/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// ---------------------------------------------------------------------------
// Potential for External Force Fields
// ---------------------------------------------------------------------------

template<typename realT>
ImplicitMembranePotential<realT>
read_implicit_membrane_potential(const toml::Table& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_implicit_membrane_potential(), 0);
    using real_type = realT;
    const auto location = "[forcefield.external] for ImplicitMembrane";

    const auto thickness = get_toml_value<real_type>(
        external, "thickness", location);
    const auto magnitude = get_toml_value<real_type>(
        external, "interaction_magnitude", location);
    const auto bend = get_toml_value<real_type>(external, "bend", location);
    MJOLNIR_LOG_INFO("thickness = ", thickness);
    MJOLNIR_LOG_INFO("magnitude = ", magnitude);
    MJOLNIR_LOG_INFO("bend      = ", bend     );

    const auto& ps = get_toml_value<toml::Array>(external, "parameters",
            location);
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    const auto parameters =
        "element of [[forcefield.external.parameters]] for ImplicitMembrane";

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index", parameters);
        const auto h   = get_toml_value<real_type>(tab, "hydrophobicity",
                                                   parameters);
        MJOLNIR_LOG_INFO("idx = ", idx, ", h = ", h);

        if(params.size() <= idx)
        {
            params.resize(idx+1, real_type(0.0));
        }
        params.at(idx) = h;
    }

    return ImplicitMembranePotential<realT>(
            thickness, magnitude, bend, std::move(params));
}

template<typename realT>
LennardJonesWallPotential<realT>
read_lennard_jones_wall_potential(const toml::Table& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_lennard_jones_wall_potential(), 0);
    using real_type = realT;
    const auto location = "[forcefield.external] for Lennard-Jones Wall";

    const auto& ps = get_toml_value<toml::Array>(external, "parameters", location);
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    const auto parameters =
        "element of [[forcefield.external.parameters]] for LennardJonesWall";

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index", parameters);
        const auto s = get_toml_value<real_type>(tab, {"sigma"_s,   u8"σ"_s}, parameters);
        const auto e = get_toml_value<real_type>(tab, {"epsilon"_s, u8"ε"_s}, parameters);
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", s, ", epsilon = ", e);

        if(params.size() <= idx)
        {
            params.resize(idx+1, std::make_pair(real_type(0.0), real_type(0.0)));
        }
        params.at(idx) = std::make_pair(s, e);
    }
    return LennardJonesWallPotential<realT>(std::move(params));
}

template<typename realT>
ExcludedVolumeWallPotential<realT>
read_excluded_volume_wall_potential(const toml::Table& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_excluded_volume_wall_potential(), 0);
    using real_type = realT;
    const auto location = "[forcefield.external] for Excluded Volume Wall";

    const real_type eps = get_toml_value<real_type>(external, "epsilon", location);
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& ps = get_toml_value<toml::Array>(external, "parameters", location);
    MJOLNIR_LOG_INFO("number of parameters = ", ps.size());

    const auto parameters =
        "element of [[forcefield.external.parameters]] for ExcludedVolumeWall";

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index", parameters);
        const auto s = get_toml_value<real_type>(tab, {"sigma"_s, u8"σ"_s}, parameters);
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", s);

        if(params.size() <= idx)
        {
            params.resize(idx+1, real_type(0.0));
        }
        params.at(idx) = s;
    }
    return ExcludedVolumeWallPotential<real_type>(eps, std::move(params));
}

} // mjolnir
#endif // MJOLNIR_READ_EXTERNAL_POTENTIAL_HPP
