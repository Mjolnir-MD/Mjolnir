#ifndef MJOLNIR_READ_POTENTIAL
#define MJOLNIR_READ_POTENTIAL
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/potential/local/Go1012ContactPotential.hpp>
#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/potential/external/ImplicitMembranePotential.hpp>
#include <mjolnir/potential/external/LennardJonesWallPotential.hpp>
#include <mjolnir/potential/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/util/get_toml_value.hpp>

namespace mjolnir
{

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, HarmonicPotential<traitsT>>>
read_harmonic_potential(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_harmonic_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, HarmonicPotential<traitsT>>;

    const auto& params = get_toml_value<toml::Array>(local, "parameters",
        "[forcefield.local] for Harmonic potential");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        using real_type = typename traitsT::real_type;
        MJOLNIR_SCOPE(for each parameters, 1);

        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = get_toml_value<indices_t>(parameter, "indices",
                "element of [[parameters]] in harmonic potential");
        auto r0 = get_toml_value<real_type>(parameter, "v0",
                "element of [[parameters]] in harmonic potential");
        auto k  = get_toml_value<real_type>(parameter, "k",
                "element of [[parameters]] in harmonic potential");
        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("r0   = ", r0);
        MJOLNIR_LOG_INFO("k    = ", k );

        retval.emplace_back(indices, HarmonicPotential<traitsT>(k, r0));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<
    std::pair<std::array<std::size_t, N>, Go1012ContactPotential<traitsT>>
    >
read_go1012_contact_potential(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_go1012_contact_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");
    if(N != 2)
    {
        MJOLNIR_LOG_WARN("Go contact potential is normally used as "
                         "a 2-body interaction");
    }

    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, Go1012ContactPotential<traitsT>>;

    const auto& params = get_toml_value<toml::Array>(local, "parameters",
            "[forcefield.local] for Go-10-12 potential");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type = typename traitsT::real_type;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = get_toml_value<indices_t>(parameter, "indices",
                "element of [[parameters]] in Go-10-12 potential");
        auto r0 = get_toml_value<real_type>(parameter, "v0",
                "element of [[parameters]] in Go-10-12 potential");
        auto k  = get_toml_value<real_type>(parameter, "k",
                "element of [[parameters]] in Go-10-12 potential");

        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("r0   = ", r0);
        MJOLNIR_LOG_INFO("k    = ", k );

        retval.emplace_back(indices, Go1012ContactPotential<traitsT>(k, r0));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, GaussianPotential<traitsT>>>
read_gaussian_potential(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_gaussian_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, GaussianPotential<traitsT>>;

    const auto& params = get_toml_value<toml::Array>(local, "parameters",
        "[forcefield.local] for Gaussian potential");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        using real_type = typename traitsT::real_type;
        MJOLNIR_SCOPE(for each parameters, 1);

        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = get_toml_value<indices_t>(parameter, "indices",
                "element of [[parameters]] in Gaussian potential");
        auto v0      = get_toml_value<real_type>(parameter, "v0",
                "element of [[parameters]] in Gaussian potential");
        auto epsilon = get_toml_value<real_type>(parameter, "epsilon",
                "element of [[parameters]] in Gaussian potential");
        auto w       = get_toml_value<real_type>(parameter, "w",
                "element of [[parameters]] in Gaussian potential");

        MJOLNIR_LOG_INFO("idxs    = ", indices);
        MJOLNIR_LOG_INFO("r0      = ", v0);
        MJOLNIR_LOG_INFO("epsilon = ", epsilon);
        MJOLNIR_LOG_INFO("w       = ", w );

        retval.emplace_back(
                indices, GaussianPotential<traitsT>(epsilon, w, v0));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<
    std::pair<std::array<std::size_t, N>, FlexibleLocalAnglePotential<traitsT>>>
read_flexible_local_angle_potential(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_flexible_local_angle_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");
    if(N != 3)
    {
        MJOLNIR_LOG_WARN("FLP angle potential is normally used as "
                         "a bond-angle interaction");
    }

    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, FlexibleLocalAnglePotential<traitsT>>;

    const auto& params = get_toml_value<toml::Array>(local, "parameters",
        "[forcefield.local] for FlexibleLocalAngle potential");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type  = typename traitsT::real_type;
        using table_type = std::array<real_type, 10>;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = get_toml_value<indices_t>(parameter, "indices",
                "element of [[parameters]] in FlexibleLocal potential");
        auto k     = get_toml_value<real_type>(parameter, "k",
                "element of [[parameters]] in FlexibleLocal potential");
        auto term1 = get_toml_value<table_type>(parameter, "term1",
                "element of [[parameters]] in FlexibleLocal potential");
        auto term2 = get_toml_value<table_type>(parameter, "term2",
                "element of [[parameters]] in FlexibleLocal potential");

        MJOLNIR_LOG_INFO("idxs  = ", indices);
        MJOLNIR_LOG_INFO("k     = ", k);
        MJOLNIR_LOG_INFO("term1 = ", term1);
        MJOLNIR_LOG_INFO("term2 = ", term2);

        retval.emplace_back(
                indices, FlexibleLocalAnglePotential<traitsT>(k, term1, term2));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>,
                      ClementiDihedralPotential<traitsT>>>
read_clementi_dihedral_potential(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_clementi_dihedral_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");
    if(N != 4)
    {
        MJOLNIR_LOG_WARN("Clementi-Go dihedral angle potential is normally used "
                         "as a dihedral-angle interaction");
    }

    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, ClementiDihedralPotential<traitsT>>;

    const auto& params = get_toml_value<toml::Array>(local, "parameters",
        "[forcefield.local] for ClementiDihedral potential");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type  = typename traitsT::real_type;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = get_toml_value<indices_t>(parameter, "indices",
                "element of [[parameters]] in ClementiDihedral potential");
        auto v0 = get_toml_value<real_type>(parameter, "v0",
                "element of [[parameters]] in ClementiDihedral potential");
        auto k1 = get_toml_value<real_type>(parameter, "k1",
                "element of [[parameters]] in ClementiDihedral potential");
        auto k3 = get_toml_value<real_type>(parameter, "k3",
                "element of [[parameters]] in ClementiDihedral potential");

        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("phi0 = ", v0);
        MJOLNIR_LOG_INFO("k1   = ", k1);
        MJOLNIR_LOG_INFO("k3   = ", k3);

        retval.emplace_back(
                indices, ClementiDihedralPotential<traitsT>(k1, k3, v0));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>,
                      FlexibleLocalDihedralPotential<traitsT>>>
read_flexible_local_dihedral_potential(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_flexible_local_dihedral_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");
    if(N != 4)
    {
        MJOLNIR_LOG_WARN("FLP angle potential is normally used as "
                         "a bond-angle interaction");
    }

    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, FlexibleLocalDihedralPotential<traitsT>>;

    const auto& params = get_toml_value<toml::Array>(local, "parameters",
        "[forcefield.local] for FlexibleLocalDihedral potential");
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type  = typename traitsT::real_type;
        using table_type = std::array<real_type, 7>;

        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = get_toml_value<indices_t>(parameter, "indices",
            "element of [[parameters]] in FlexibleLocalDihedral potential");
        auto k    = get_toml_value<real_type>(parameter, "k",
            "element of [[parameters]] in FlexibleLocalDihedral potential");
        auto term = get_toml_value<table_type>(parameter, "term",
            "element of [[parameters]] in FlexibleLocalDihedral potential");

        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("k    = ", k);
        MJOLNIR_LOG_INFO("term = ", term);

       retval.emplace_back(
                indices, FlexibleLocalDihedralPotential<traitsT>(k, term));
    }
    return retval;
}

// ----------------------------------------------------------------------------
// global potential
// ----------------------------------------------------------------------------

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
        MJOLNIR_SCOPE(for each parameters, 1);

        const auto& tab = param.cast<toml::value_t::Table>();
        const auto  idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.global.parameters]] for Excluded Volume");
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }

        const auto sgm = get_toml_value<real_type>(tab, "sigma",
            "element of [[forcefield.global.parameters]] for Excluded Volume");
        MJOLNIR_LOG_INFO("idx = ", idx);
        MJOLNIR_LOG_INFO("sgm = ", sgm);
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
        MJOLNIR_SCOPE(for each parameters, 1);
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
        MJOLNIR_LOG_INFO("idx = ", idx);
        MJOLNIR_LOG_INFO("sgm = ", sigma);
        MJOLNIR_LOG_INFO("eps = ", epsilon);

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
        MJOLNIR_SCOPE(for each parameters, 1);
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.global.parameters]] for Debye-Huckel");
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        const auto charge = get_toml_value<real_type>(tab, "charge",
            "element of [[forcefield.global.parameters]] for Debye-Huckel");
        MJOLNIR_LOG_INFO("idx    = ", idx);
        MJOLNIR_LOG_INFO("charge = ", charge);
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
        MJOLNIR_SCOPE(for each parameters, 1);
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
        MJOLNIR_LOG_INFO("idx = ", idx);
        MJOLNIR_LOG_INFO("h   = ", h);
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
        MJOLNIR_SCOPE(for each parameters, 1);
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
        MJOLNIR_LOG_INFO("idx     = ", idx);
        MJOLNIR_LOG_INFO("sigma   = ", s);
        MJOLNIR_LOG_INFO("epsilon = ", e);
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
        MJOLNIR_SCOPE(for each parameters, 1);
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = get_toml_value<std::size_t>(tab, "index",
            "element of [[forcefield.external.parameters]] for ExcludedVolumeWall");
        if(params.size() <= idx)
        {
            params.resize(idx+1, real_type(0.0));
        }

        const auto s = get_toml_value<real_type>(tab, "sigma",
            "element of [[forcefield.external.parameters]] for ExcludedVolumeWall");
        MJOLNIR_LOG_INFO("idx   = ", idx);
        MJOLNIR_LOG_INFO("sigma = ", s);
        params.at(idx) = s;
    }
    return ExcludedVolumeWallPotential<traitsT>(eps, std::move(params));
}


} // mjolnir
#endif // MJOLNIR_READ_POTENTIAL
