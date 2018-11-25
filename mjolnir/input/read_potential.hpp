#ifndef MJOLNIR_READ_POTENTIAL
#define MJOLNIR_READ_POTENTIAL
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/potential/local/Go1012ContactPotential.hpp>
#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>
#include <mjolnir/potential/local/AngularGaussianPotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/potential/local/SumLocalPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
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

// this reads...
// parameters = [
// {indices = [1, 2], k = 10.0, v0 = 1.0} <- a portion of this table, k and v0.
// ]
template<typename realT>
HarmonicPotential<realT> read_harmonic_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto location = "element of [[parameters]] in harmonic potential";

    const auto k  = get_toml_value<real_type>(param, "k",  location);
    const auto v0 = get_toml_value<real_type>(param, "v0", location);

    MJOLNIR_LOG_INFO("HarmonicPotential = {v0 = ", v0, ", k = ", k, '}');
    return HarmonicPotential<realT>(k, v0);
}

template<typename realT>
Go1012ContactPotential<realT>
read_go1012_contact_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto location = "element of [[parameters]] in Go-10-12 potential";

    const auto k  = get_toml_value<real_type>(param, "k",  location);
    const auto v0 = get_toml_value<real_type>(param, "v0", location);

    MJOLNIR_LOG_INFO("GoContactPotential = {v0 = ", v0, ", k = ", k, '}');
    return Go1012ContactPotential<realT>(k, v0);
}

template<typename realT>
GaussianPotential<realT> read_gaussian_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto location = "element of [[parameters]] in Gaussian potential";

    const auto v0    = get_toml_value<real_type>(param, "v0",    location);
    const auto k     = get_toml_value<real_type>(param, "k",     location);
    const auto sigma =
        get_toml_value<real_type>(param, {"sigma"_s, u8"σ"_s}, location);

    MJOLNIR_LOG_INFO("GaussianPotential = {v0 = ", v0, ", k = ", k,
                     ", sigma = ", sigma, '}');
    return GaussianPotential<realT>(k, sigma, v0);
}

template<typename realT>
AngularGaussianPotential<realT>
read_angular_gaussian_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto location =
        "element of [[parameters]] in Angular Gaussian potential";

    const auto v0    = get_toml_value<real_type>(param, "v0",    location);
    const auto k     = get_toml_value<real_type>(param, "k",     location);
    const auto sigma =
        get_toml_value<real_type>(param, {"sigma"_s, u8"σ"_s}, location);

    MJOLNIR_LOG_INFO("AngularGaussianPotential = {v0 = ", v0, ", k = ", k,
                     ", sigma = ", sigma, '}');
    return AngularGaussianPotential<realT>(k, sigma, v0);
}

template<typename realT>
FlexibleLocalAnglePotential<realT>
read_flexible_local_angle_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    using term_type = std::array<real_type, 10>;
    const auto location =
        "element of [[parameters]] in FlexibleLocal potential";

    const auto k     = get_toml_value<real_type>(param, "k",     location);
    const auto term1 = get_toml_value<term_type>(param, "term1", location);
    const auto term2 = get_toml_value<term_type>(param, "term2", location);

    MJOLNIR_LOG_INFO("FlexibleLocalAngle = {k = ", k,
                     ", term1 = ", term1, ", term2 = ", term2, '}');
    return FlexibleLocalAnglePotential<realT>(k, term1, term2);
}

template<typename realT>
ClementiDihedralPotential<realT>
read_clementi_dihedral_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto location =
        "element of [[parameters]] in ClementiDihedral potential";

    const auto v0 = get_toml_value<real_type>(param, "v0", location);
    const auto k1 = get_toml_value<real_type>(param, "k1", location);
    const auto k3 = get_toml_value<real_type>(param, "k3", location);

    MJOLNIR_LOG_INFO("ClementiDihedral = {v0 = ", v0,
                     ", k1 = ", k1, ", k3 = ", k3, '}');
    return ClementiDihedralPotential<realT>(k1, k3, v0);
}

template<typename realT>
FlexibleLocalDihedralPotential<realT>
read_flexible_local_dihedral_potential(const toml::Table& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    using term_type = std::array<real_type, 7>;
    const auto location =
        "element of [[parameters]] in FlexibleLocalDihedral potential";

    auto k    = get_toml_value<real_type>(param, "k",    location);
    auto term = get_toml_value<term_type>(param, "term", location);

    MJOLNIR_LOG_INFO("FlexibleLocalDihedral = {k = ", k, ", term = ", term, '}');
    return FlexibleLocalDihedralPotential<realT>(k, term);
}

// ----------------------------------------------------------------------------
// utility function to read local potentials
// ----------------------------------------------------------------------------

// potential_class -> reading_function adapter to implement a function
// `read_local_potential` that reads indices and call read_*_potential.
// local potential requires not only the function parameter but also the indices
// of particles. to split them, mjolnir has both `read_local_potential` and
// `read_(harmonic|gaussian|...)_potential`.
namespace detail
{
template<typename potentialT> struct read_local_potential_impl;

template<typename realT>
struct read_local_potential_impl<HarmonicPotential<realT>>
{
    static HarmonicPotential<realT> invoke(const toml::Table& param)
    {
        return read_harmonic_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<Go1012ContactPotential<realT>>
{
    static Go1012ContactPotential<realT> invoke(const toml::Table& param)
    {
        return read_go1012_contact_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<GaussianPotential<realT>>
{
    static GaussianPotential<realT> invoke(const toml::Table& param)
    {
        return read_gaussian_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<AngularGaussianPotential<realT>>
{
    static AngularGaussianPotential<realT> invoke(const toml::Table& param)
    {
        return read_angular_gaussian_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<FlexibleLocalAnglePotential<realT>>
{
    static FlexibleLocalAnglePotential<realT> invoke(const toml::Table& param)
    {
        return read_flexible_local_angle_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<ClementiDihedralPotential<realT>>
{
    static ClementiDihedralPotential<realT> invoke(const toml::Table& param)
    {
        return read_clementi_dihedral_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<FlexibleLocalDihedralPotential<realT>>
{
    static FlexibleLocalDihedralPotential<realT> invoke(const toml::Table& param)
    {
        return read_flexible_local_dihedral_potential<realT>(param);
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
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

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

template<std::size_t N, typename realT,
         typename potentialT1, typename potentialT2>
std::vector<std::pair<std::array<std::size_t, N>,
            SumLocalPotential<realT, potentialT1, potentialT2>>>
read_local_potentials(const toml::Table& local,
                      const std::string& p1, const std::string& p2)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_local_potentials(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t                = std::array<std::size_t, N>;
    using indices_potential_pair_t = std::pair<indices_t,
        SumLocalPotential<realT, potentialT1, potentialT2>>;

    const auto& params = get_toml_value<toml::Array>(
            local, "parameters", "[[forcefield.local]]");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

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
            SumLocalPotential<realT, potentialT1, potentialT2>(
                detail::read_local_potential_impl<potentialT1>::invoke(pot1),
                detail::read_local_potential_impl<potentialT2>::invoke(pot2)
            ));
    }
    return retval;
}


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
#endif // MJOLNIR_READ_POTENTIAL
