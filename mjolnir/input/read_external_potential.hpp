#ifndef MJOLNIR_READ_EXTERNAL_POTENTIAL_HPP
#define MJOLNIR_READ_EXTERNAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/external/ImplicitMembranePotential.hpp>
#include <mjolnir/potential/external/LennardJonesWallPotential.hpp>
#include <mjolnir/potential/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/core/Topology.hpp>
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
read_implicit_membrane_potential(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_implicit_membrane_potential(), 0);
    using real_type = realT;

    const auto thickness = toml::find<real_type>(external, "thickness");
    const auto magnitude = toml::find<real_type>(external, "interaction_magnitude");
    const auto bend      = toml::find<real_type>(external, "bend");

    MJOLNIR_LOG_INFO("thickness = ", thickness);
    MJOLNIR_LOG_INFO("magnitude = ", magnitude);
    MJOLNIR_LOG_INFO("bend      = ", bend     );

    const auto& ps = toml::find<toml::array>(external, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = toml::find<std::size_t>(ps, "index");
        const auto h   = toml::find<real_type  >(ps, "hydrophobicity");
        if(params.size() <= idx) {params.resize(idx+1, real_type(0.0));}
        params.at(idx) = h;

        MJOLNIR_LOG_INFO("idx = ", idx, ", hydrophobicity = ", h);
    }
    return ImplicitMembranePotential<realT>(
            thickness, magnitude, bend, std::move(params));
}

template<typename realT>
LennardJonesWallPotential<realT>
read_lennard_jones_wall_potential(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_lennard_jones_wall_potential(), 0);
    using real_type = realT;

    const auto& ps = toml::find<toml::array>(external, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = toml::find<std::size_t>(ps, "index");
        const auto s   = toml::expect<real_type>(ps, u8"σ").or_other(
                         toml::expect<real_type>(ps, "sigma")).unwrap();
        const auto e   = toml::expect<real_type>(ps, u8"ε").or_other(
                         toml::expect<real_type>(ps, "epsilon")).unwrap();
        if(params.size() <= idx)
        {
            params.resize(idx+1, std::make_pair(real_type(0), real_type(0)));
        }
        params.at(idx) = std::make_pair(s, e);

        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", s, ", epsilon = ", e);
    }
    return LennardJonesWallPotential<realT>(std::move(params));
}

template<typename realT>
ExcludedVolumeWallPotential<realT>
read_excluded_volume_wall_potential(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_excluded_volume_wall_potential(), 0);
    using real_type = realT;

    const auto eps = toml::expect<real_type>(external, u8"ε").or_other(
                     toml::expect<real_type>(external, "epsilon")).unwrap();
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& ps = toml::find<toml::array>(external, "parameters");
    MJOLNIR_LOG_INFO("number of parameters = ", ps.size());

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = toml::find<std::size_t>(ps, "index");
        const auto s   = toml::expect<real_type>(ps, u8"σ").or_other(
                         toml::expect<real_type>(ps, "sigma")).unwrap();
        if(params.size() <= idx) {params.resize(idx+1, real_type(0));}
        params.at(idx) = s;

        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", s);
    }
    return ExcludedVolumeWallPotential<real_type>(eps, std::move(params));
}

} // mjolnir
#endif // MJOLNIR_READ_EXTERNAL_POTENTIAL_HPP
