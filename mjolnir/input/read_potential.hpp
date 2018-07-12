#ifndef MJOLNIR_READ_POTENTIAL
#define MJOLNIR_READ_POTENTIAL
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/potential/local/Go1012ContactPotential.hpp>
#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/potential/local/SumLocalPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/potential/external/ImplicitMembranePotential.hpp>
#include <mjolnir/potential/external/LennardJonesWallPotential.hpp>
#include <mjolnir/potential/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// this reads...
// parameters = [
// {indices = [1, 2], k = 10.0, v0 = 1.0} <- a portion of this table, k and v0.
// ]
template<typename traitsT>
HarmonicPotential<traitsT> read_harmonic_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;
    const auto location = "element of [[parameters]] in harmonic potential";

    const auto k  = get_toml_value<real_type>(param, "k",  location);
    const auto v0 = get_toml_value<real_type>(param, "v0", location);

    MJOLNIR_LOG_INFO("HarmonicPotential = {v0 = ", v0, ", k = ", k, '}');
    return HarmonicPotential<traitsT>(k, v0);
}

template<typename traitsT>
Go1012ContactPotential<traitsT>
read_go1012_contact_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;
    const auto location = "element of [[parameters]] in Go-10-12 potential";

    const auto k  = get_toml_value<real_type>(param, "k",  location);
    const auto v0 = get_toml_value<real_type>(param, "v0", location);

    MJOLNIR_LOG_INFO("GoContactPotential = {v0 = ", v0, ", k = ", k, '}');
    return Go1012ContactPotential<traitsT>(k, v0);
}

template<typename traitsT>
GaussianPotential<traitsT> read_gaussian_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;
    const auto location = "element of [[parameters]] in Gaussian potential";

    const auto v0      = get_toml_value<real_type>(param, "v0",      location);
    const auto epsilon = get_toml_value<real_type>(param, "epsilon", location);
    const auto w       = get_toml_value<real_type>(param, "w",       location);

    MJOLNIR_LOG_INFO("GaussianPotential = {v0 = ", v0,
                     ", epsilon = ", epsilon, ", w = ", w, '}');
    return GaussianPotential<traitsT>(epsilon, w, v0);
}

template<typename traitsT>
FlexibleLocalAnglePotential<traitsT>
read_flexible_local_angle_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;
    using term_type = std::array<real_type, 10>;
    const auto location =
        "element of [[parameters]] in FlexibleLocal potential";

    const auto k     = get_toml_value<real_type>(param, "k",     location);
    const auto term1 = get_toml_value<term_type>(param, "term1", location);
    const auto term2 = get_toml_value<term_type>(param, "term2", location);

    MJOLNIR_LOG_INFO("FlexibleLocalAngle = {k = ", k,
                     ", term1 = ", term1, ", term2 = ", term2, '}');
    return FlexibleLocalAnglePotential<traitsT>(k, term1, term2);
}

template<typename traitsT>
ClementiDihedralPotential<traitsT>
read_clementi_dihedral_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type  = typename traitsT::real_type;
    const auto location =
        "element of [[parameters]] in ClementiDihedral potential";

    const auto v0 = get_toml_value<real_type>(param, "v0", location);
    const auto k1 = get_toml_value<real_type>(param, "k1", location);
    const auto k3 = get_toml_value<real_type>(param, "k3", location);

    MJOLNIR_LOG_INFO("ClementiDihedral = {v0 = ", v0,
                     ", k1 = ", k1, ", k3 = ", k3, '}');
    return ClementiDihedralPotential<traitsT>(k1, k3, v0);
}

template<typename traitsT>
FlexibleLocalDihedralPotential<traitsT>
read_flexible_local_dihedral_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;
    using term_type = std::array<real_type, 7>;
    const auto location =
        "element of [[parameters]] in FlexibleLocalDihedral potential";

    auto k    = get_toml_value<real_type>(param, "k",    location);
    auto term = get_toml_value<term_type>(param, "term", location);

    MJOLNIR_LOG_INFO("FlexibleLocalDihedral = {k = ", k, ", term = ", term, '}');
    return FlexibleLocalDihedralPotential<traitsT>(k, term);
}

// ----------------------------------------------------------------------------
// utility function to read local potentials
// ----------------------------------------------------------------------------

// potential_class -> reading_function adapter
// to enable partial specialization of function template, mjolnir uses
// static function inside of struct. struct can be partially specialized.
namespace detail
{
template<typename potentialT> struct read_local_potential_impl;

template<typename traitsT>
struct read_local_potential_impl<HarmonicPotential<traitsT>>
{
    static HarmonicPotential<traitsT> invoke(const toml::Table& param)
    {
        return read_harmonic_potential<traitsT>(param);
    }
};
template<typename traitsT>
struct read_local_potential_impl<Go1012ContactPotential<traitsT>>
{
    static Go1012ContactPotential<traitsT> invoke(const toml::Table& param)
    {
        return read_go1012_contact_potential<traitsT>(param);
    }
};
template<typename traitsT>
struct read_local_potential_impl<GaussianPotential<traitsT>>
{
    static GaussianPotential<traitsT> invoke(const toml::Table& param)
    {
        return read_gaussian_potential<traitsT>(param);
    }
};
template<typename traitsT>
struct read_local_potential_impl<FlexibleLocalAnglePotential<traitsT>>
{
    static FlexibleLocalAnglePotential<traitsT> invoke(const toml::Table& param)
    {
        return read_flexible_local_angle_potential<traitsT>(param);
    }
};
template<typename traitsT>
struct read_local_potential_impl<ClementiDihedralPotential<traitsT>>
{
    static ClementiDihedralPotential<traitsT> invoke(const toml::Table& param)
    {
        return read_clementi_dihedral_potential<traitsT>(param);
    }
};
template<typename traitsT>
struct read_local_potential_impl<FlexibleLocalDihedralPotential<traitsT>>
{
    static FlexibleLocalDihedralPotential<traitsT> invoke(const toml::Table& param)
    {
        return read_flexible_local_dihedral_potential<traitsT>(param);
    }
};
} // namespace detail

// this function reads particle indices on which the potential will be applied
// and returns pairs of [indices, potential parameters].
template<std::size_t N, typename potentialT>
std::vector<std::pair<std::array<std::size_t, N>, potentialT>>
read_local_potential(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_local_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t                = std::array<std::size_t, N>;
    using indices_potential_pair_t = std::pair<indices_t, potentialT>;

    const auto& params = get_toml_value<toml::Array>(
            local, "parameters", "[[forcefield.local]]");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();

        const auto indices = get_toml_value<indices_t>(parameter, "indices",
                "element of [[forcefields.local.parameters]]");
        MJOLNIR_LOG_INFO_NO_LF("idxs = ", indices, ", ");

        retval.emplace_back(indices,
            detail::read_local_potential_impl<potentialT>::invoke(parameter));
    }
    return retval;
}

template<std::size_t N, typename traitsT,
         typename potentialT1, typename potentialT2>
std::vector<std::pair<std::array<std::size_t, N>,
            SumLocalPotential<traitsT, potentialT1, potentialT2>>>
read_local_potentials(const toml::Table& local,
                      const std::string& p1, const std::string& p2)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_local_potentials(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t                = std::array<std::size_t, N>;
    using indices_potential_pair_t = std::pair<indices_t,
        SumLocalPotential<traitsT, potentialT1, potentialT2>>;

    const auto& params = get_toml_value<toml::Array>(
            local, "parameters", "[[forcefield.local]]");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();

        const auto indices = get_toml_value<indices_t>(parameter, "indices",
                "element of [[forcefields.local.parameters]]");
        MJOLNIR_LOG_INFO("idxs = ", indices);

        const auto& pot1 =
            get_toml_value<toml::Table>(parameter, p1, "[[forcefield.local]]");
        const auto& pot2 =
            get_toml_value<toml::Table>(parameter, p2, "[[forcefield.local]]");

        retval.emplace_back(indices,
            SumLocalPotential<traitsT, potentialT1, potentialT2>(
                detail::read_local_potential_impl<potentialT1>::invoke(pot1),
                detail::read_local_potential_impl<potentialT2>::invoke(pot2)
            ));
    }
    return retval;
}


// ============================================================================
// global potential
// ============================================================================

template<typename traitsT, typename ignoreT>
ExcludedVolumePotential<traitsT, ignoreT>
read_excluded_volume_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_excluded_volume_potential(), 0);
    typedef typename traitsT::real_type real_type;

    const auto& ignored_connections = get_toml_value<toml::Table>(
        global, "ignored_connections", "[forcefield.global] for ExcludedVolume");
    std::map<std::string, std::size_t> connections;
    for(const auto connection : ignored_connections)
    {
        connections[connection.first] =
            toml::get<std::size_t>(connection.second);
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
                         " within ", connections.at(connection.first), " will ",
                         "be ignored");
    }

    const real_type eps = get_toml_value<real_type>(
        global, "epsilon", "[forcefield.global] for ExcludedVolume");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& ps = get_toml_value<toml::Array>(global, "parameters",
            "[forcefield.global] for ExcludedVolume potential");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());

    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto  idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.global.parameters]] for Excluded Volume");
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }

        const auto sgm = get_toml_value<real_type>(tab, "sigma",
            "element of [[forcefield.global.parameters]] for Excluded Volume");
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sgm);
        params.at(idx) = sgm;
    }

    return ExcludedVolumePotential<traitsT, ignoreT>(
            eps, std::move(params), connections);
}

template<typename traitsT, typename ignoreT>
LennardJonesPotential<traitsT, ignoreT>
read_lennard_jones_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_lennard_jones_potential(), 0);
    typedef typename traitsT::real_type real_type;

    const auto& ignored_connections = get_toml_value<toml::Table>(
        global, "ignored_connections", "[forcefield.global] for ExcludedVolume");
    std::map<std::string, std::size_t> connections;
    for(const auto connection : ignored_connections)
    {
        connections[connection.first] =
            toml::get<std::size_t>(connection.second);
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
                         " within ", connections.at(connection.first), " will ",
                         "be ignored");
    }

    const auto& ps = get_toml_value<toml::Array>(global, "parameters",
        "[forcefield.global] for Lennard-Jones potential");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.global.parameters]] for Lennard-Jones");
        if(params.size() <= idx)
        {
            const std::pair<real_type, real_type> dummy{0., 0.};
            params.resize(idx+1, dummy);
        }
        const auto sigma   = get_toml_value<real_type>(tab, "sigma",
            "element of [[forcefield.global.parameters]] for Lennard-Jones");
        const auto epsilon = get_toml_value<real_type>(tab, "epsilon",
            "element of [[forcefield.global.parameters]] for Lennard-Jones");
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma,
                         ", epsilon = ", epsilon);

        params.at(idx) = std::make_pair(sigma, epsilon);
    }

    return LennardJonesPotential<traitsT, ignoreT>(
            std::move(params), connections);
}

template<typename traitsT, typename ignoreT>
DebyeHuckelPotential<traitsT, ignoreT>
read_debye_huckel_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_debye_huckel_potential(), 0);
    typedef typename traitsT::real_type real_type;

    const auto& ignored_connections = get_toml_value<toml::Table>(
        global, "ignored_connections", "[forcefield.global] for ExcludedVolume");
    std::map<std::string, std::size_t> connections;
    for(const auto connection : ignored_connections)
    {
        connections[connection.first] =
            toml::get<std::size_t>(connection.second);
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
                         " within ", connections.at(connection.first), " will ",
                         "be ignored");
    }

    const auto& ps = get_toml_value<toml::Array>(global, "parameters",
        "[forcefield.global] for Debye-Huckel");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.global.parameters]] for Debye-Huckel");
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        const auto charge = get_toml_value<real_type>(tab, "charge",
            "element of [[forcefield.global.parameters]] for Debye-Huckel");
        MJOLNIR_LOG_INFO("idx    = ", idx, ", charge = ", charge);
        params.at(idx) = charge;
    }

    return DebyeHuckelPotential<traitsT, ignoreT>(
            std::move(params), connections);
}

// ---------------------------------------------------------------------------
// Potential for External Force Fields
// ---------------------------------------------------------------------------

template<typename traitsT>
ImplicitMembranePotential<traitsT>
read_implicit_membrane_potential(const toml::Table& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_implicit_membrane_potential(), 0);
    typedef typename traitsT::real_type real_type;

    const auto thickness = get_toml_value<real_type>(external,
        "thickness", "[forcefield.external] for ImplicitMembrane");
    const auto magnitude = get_toml_value<real_type>(external,
        "interaction_magnitude", "[forcefield.external] for ImplicitMembrane");
    const auto bend = get_toml_value<real_type>(external,
        "bend", "[forcefield.external] for ImplicitMembrane");
    MJOLNIR_LOG_INFO("thickness = ", thickness);
    MJOLNIR_LOG_INFO("magnitude = ", magnitude);
    MJOLNIR_LOG_INFO("bend      = ", bend     );

    const auto& ps = get_toml_value<toml::Array>(external, "parameters",
            "[forcefield.external] for ImplicitMembrane");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.external.parameters]] for ImplicitMembrane");

        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        const auto h = get_toml_value<real_type>(tab, "hydrophobicity",
            "element of [[forcefield.external.parameters]] for ImplicitMembrane"
            );
        MJOLNIR_LOG_INFO("idx = ", idx, ", h = ", h);
        params.at(idx) = h;
    }

    return ImplicitMembranePotential<traitsT>(
            thickness, magnitude, bend, std::move(params));
}

template<typename traitsT>
LennardJonesWallPotential<traitsT>
read_lennard_jones_wall_potential(const toml::Table& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_lennard_jones_wall_potential(), 0);
    typedef typename traitsT::real_type real_type;

    const auto& ps = get_toml_value<toml::Array>(external, "parameters",
            "[forcefield.external] for Lennard-Jones Wall");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.external.parameters]] for LennardJonesWall");
        if(params.size() <= idx)
        {
            params.resize(idx+1, std::make_pair(real_type(0.0), real_type(0.0)));
        }

        const auto s = get_toml_value<real_type>(tab, "sigma",
            "element of [[forcefield.external.parameters]] for LennardJonesWall");
        const auto e = get_toml_value<real_type>(tab, "epsilon",
            "element of [[forcefield.external.parameters]] for LennardJonesWall");
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", s, ", epsilon = ", e);
        params.at(idx) = std::make_pair(s, e);
    }
    return LennardJonesWallPotential<traitsT>(std::move(params));
}

template<typename traitsT>
ExcludedVolumeWallPotential<traitsT>
read_excluded_volume_wall_potential(const toml::Table& external)
{
    typedef typename traitsT::real_type real_type;
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_excluded_volume_wall_potential(), 0);

    const auto& ps = get_toml_value<toml::Array>(external, "parameters",
            "[forcefield.external] for Excluded Volume Wall");
    MJOLNIR_LOG_INFO("number of parameters = ", ps.size());

    const real_type eps = get_toml_value<real_type>(
        external, "epsilon", "[forcefield.external] for ExcludedVolume");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.external.parameters]] for ExcludedVolumeWall");
        if(params.size() <= idx)
        {
            params.resize(idx+1, real_type(0.0));
        }

        const auto s = get_toml_value<real_type>(tab, "sigma",
            "element of [[forcefield.external.parameters]] for ExcludedVolumeWall");
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", s);
        params.at(idx) = s;
    }
    return ExcludedVolumeWallPotential<traitsT>(eps, std::move(params));
}


} // mjolnir
#endif // MJOLNIR_READ_POTENTIAL
