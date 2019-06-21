#ifndef MJOLNIR_INPUT_READ_EXTERNAL_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_EXTERNAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/potential/external/HarmonicRestraintPotential.hpp>
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
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto thickness = toml::find<real_type>(external, "thickness");
    const auto magnitude = toml::find<real_type>(external, "interaction_magnitude");
    const auto bend      = toml::find<real_type>(external, "bend");

    MJOLNIR_LOG_INFO("thickness = ", thickness);
    MJOLNIR_LOG_INFO("magnitude = ", magnitude);
    MJOLNIR_LOG_INFO("bend      = ", bend     );

    const auto& env = external.as_table().count("env") == 1 ?
                      external.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(external, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index");
        const auto h   = find_parameter<real_type  >(param, env, "hydrophobicity");
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
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& env = external.as_table().count("env") == 1 ?
                      external.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(external, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index");
        const auto s   = find_parameter<real_type>(param, env, "sigma", u8"σ");
        const auto e   = find_parameter<real_type>(param, env, "epsilon", u8"ε");
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
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto eps = toml::expect<real_type>(external, u8"ε").or_other(
                     toml::expect<real_type>(external, "epsilon")).unwrap();
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& env = external.as_table().count("env") == 1 ?
                      external.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(external, "parameters");
    MJOLNIR_LOG_INFO("number of parameters = ", ps.size());

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index");
        const auto rad = find_parameter<real_type  >(param, env, "radius");
        if(params.size() <= idx) {params.resize(idx+1, real_type(0));}
        params.at(idx) = rad;

        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", rad);
    }
    return ExcludedVolumeWallPotential<real_type>(eps, std::move(params));
}

template<typename realT>
HarmonicRestraintPotential<realT>
read_harmonic_restraint_potential(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;
    using potential_type = HarmonicRestraintPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;

    const auto& ps = toml::find<toml::array>(external, "parameters");
    MJOLNIR_LOG_INFO("number of parameters = ", ps.size());

    const auto& env = external.as_table().count("env") == 1 ?
                      external.as_table().at("env") : toml::value{};

    std::vector<parameter_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index");
        const auto k   = find_parameter<real_type  >(param, env, "k");
        const auto v0  = find_parameter<real_type  >(param, env, "v0");
        if(params.size() <= idx)
        {
            params.resize(idx+1, parameter_type(0, 0));
        }
        params.at(idx) = parameter_type(k, v0);

        MJOLNIR_LOG_INFO("idx = ", idx, ", k = ", k, ", v0 = ", v0);
    }
    return HarmonicRestraintPotential<real_type>(std::move(params));
}

} // mjolnir
#endif // MJOLNIR_READ_EXTERNAL_POTENTIAL_HPP
