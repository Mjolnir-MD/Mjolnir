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

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for Harmonic potential"
        ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    MJOLNIR_LOG_INFO("{{{");
    for(const auto& item : params)
    {
        using real_type = typename traitsT::real_type;
        MJOLNIR_SCOPE(for each parameters, 1);

        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in harmonic potential"));
        auto r0 = toml::get<real_type>(toml_value_at(parameter, "eq",
                "element of [[parameters]] in harmonic potential"));
        auto k  = toml::get<real_type>(toml_value_at(parameter, "k",
                "element of [[parameters]] in harmonic potential"));
        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("r0   = ", r0);
        MJOLNIR_LOG_INFO("k    = ", k );

        retval.emplace_back(indices, HarmonicPotential<traitsT>(k, r0));
    }
    MJOLNIR_LOG_INFO("}}}");
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

    const auto& params = toml_value_at(local, "parameters",
            "[forcefield.local] for Go-10-12 potential"
            ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    MJOLNIR_LOG_INFO("{{{");
    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type = typename traitsT::real_type;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in Go-10-12 potential"));
        auto r0 = toml::get<real_type>(toml_value_at(parameter, "eq",
                "element of [[parameters]] in Go-10-12 potential"));
        auto k  = toml::get<real_type>(toml_value_at(parameter, "k",
                "element of [[parameters]] in Go-10-12 potential"));

        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("r0   = ", r0);
        MJOLNIR_LOG_INFO("k    = ", k );

        retval.emplace_back(indices, Go1012ContactPotential<traitsT>(k, r0));
    }
    MJOLNIR_LOG_INFO("}}}");
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

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for Gaussian potential"
        ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    MJOLNIR_LOG_INFO("{{{");
    for(const auto& item : params)
    {
        using real_type = typename traitsT::real_type;
        MJOLNIR_SCOPE(for each parameters, 1);

        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in Gaussian potential"));
        auto eq      = toml::get<real_type>(toml_value_at(parameter, "eq",
                "element of [[parameters]] in Gaussian potential"));
        auto epsilon = toml::get<real_type>(toml_value_at(parameter, "epsilon",
                "element of [[parameters]] in Gaussian potential"));
        auto w       = toml::get<real_type>(toml_value_at(parameter, "w",
                "element of [[parameters]] in Gaussian potential"));

        MJOLNIR_LOG_INFO("idxs    = ", indices);
        MJOLNIR_LOG_INFO("r0      = ", eq);
        MJOLNIR_LOG_INFO("epsilon = ", epsilon);
        MJOLNIR_LOG_INFO("w       = ", w );

        retval.emplace_back(
                indices, GaussianPotential<traitsT>(epsilon, w, eq));
    }
    MJOLNIR_LOG_INFO("}}}");
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

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for FlexibleLocalAngle potential"
        ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    MJOLNIR_LOG_INFO("{{{");
    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type  = typename traitsT::real_type;
        using table_type = std::array<real_type, 10>;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in FlexibleLocal potential"));
        auto k     = toml::get<real_type>(toml_value_at(parameter, "k",
                "element of [[parameters]] in FlexibleLocal potential"));
        auto term1 = toml::get<table_type>(toml_value_at(parameter, "term1",
                "element of [[parameters]] in FlexibleLocal potential"));
        auto term2 = toml::get<table_type>(toml_value_at(parameter, "term2",
                "element of [[parameters]] in FlexibleLocal potential"));

        MJOLNIR_LOG_INFO("idxs  = ", indices);
        MJOLNIR_LOG_INFO("k     = ", k);
        MJOLNIR_LOG_INFO("term1 = ", term1);
        MJOLNIR_LOG_INFO("term2 = ", term2);

        retval.emplace_back(
                indices, FlexibleLocalAnglePotential<traitsT>(k, term1, term2));
    }
    MJOLNIR_LOG_INFO("}}}");
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

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for ClementiDihedral potential"
        ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    MJOLNIR_LOG_INFO("{{{");
    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type  = typename traitsT::real_type;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in ClementiDihedral potential"));
        auto eq = toml::get<real_type>(toml_value_at(parameter, "eq",
                "element of [[parameters]] in ClementiDihedral potential"));
        auto k1 = toml::get<real_type>(toml_value_at(parameter, "k1",
                "element of [[parameters]] in ClementiDihedral potential"));
        auto k3 = toml::get<real_type>(toml_value_at(parameter, "k3",
                "element of [[parameters]] in ClementiDihedral potential"));

        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("phi0 = ", eq);
        MJOLNIR_LOG_INFO("k1   = ", k1);
        MJOLNIR_LOG_INFO("k3   = ", k3);

        retval.emplace_back(
                indices, ClementiDihedralPotential<traitsT>(k1, k3, eq));
    }
    MJOLNIR_LOG_INFO("}}}");
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

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for FlexibleLocalDihedral potential"
        ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(params.size(), " parameters are found");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    MJOLNIR_LOG_INFO("{{{");
    for(const auto& item : params)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        using real_type  = typename traitsT::real_type;
        using table_type = std::array<real_type, 7>;

        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
            "element of [[parameters]] in FlexibleLocalDihedral potential"));
        auto k    = toml::get<real_type>(toml_value_at(parameter, "k",
            "element of [[parameters]] in FlexibleLocalDihedral potential"));
        auto term = toml::get<table_type>(toml_value_at(parameter, "term",
            "element of [[parameters]] in FlexibleLocalDihedral potential"));

        MJOLNIR_LOG_INFO("idxs = ", indices);
        MJOLNIR_LOG_INFO("k    = ", k);
        MJOLNIR_LOG_INFO("term = ", term);

       retval.emplace_back(
                indices, FlexibleLocalDihedralPotential<traitsT>(k, term));
    }
    MJOLNIR_LOG_INFO("}}}");
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

    const std::size_t bonds = toml::get<std::size_t>(toml_value_at(
        global, "ignored_bonds", "[forcefield.global] for ExcludedVolume"));
    const std::size_t contacts = toml::get<std::size_t>(toml_value_at(
        global, "ignored_contacts", "[forcefield.global] for ExcludedVolume"));
    const real_type eps = toml::get<real_type>(toml_value_at(
        global, "epsilon", "[forcefield.global] for ExcludedVolume"));
    MJOLNIR_LOG_INFO("particles connected less than ", bonds,    " bonds are ignored");
    MJOLNIR_LOG_INFO("particles connected less than ", contacts, " contacts are ignored");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& ps = toml_value_at(global, "parameters",
            "[forcefield.global] for ExcludedVolume potential"
            ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());

    MJOLNIR_LOG_INFO("{{{");
    for(const auto& param : ps)
    {
        MJOLNIR_SCOPE(for each parameters, 1);

        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.global.parameters]] for Excluded Volume"));
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }

        const auto sgm = toml::get<real_type>(toml_value_at(tab, "sigma",
            "element of [[forcefield.global.parameters]] for Excluded Volume"));
        MJOLNIR_LOG_INFO("idx = ", idx);
        MJOLNIR_LOG_INFO("sgm = ", sgm);
        params.at(idx) = sgm;
    }
    MJOLNIR_LOG_INFO("}}}");

    return ExcludedVolumePotential<traitsT, ignoreT>(
            eps, std::move(params), bonds, contacts);
}

template<typename traitsT, typename ignoreT>
LennardJonesPotential<traitsT, ignoreT>
read_lennard_jones_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_lennard_jones_potential(), 0);
    typedef typename traitsT::real_type real_type;

    const std::size_t bonds = toml::get<std::size_t>(toml_value_at(
        global, "ignored_bonds", "[forcefield.global] for Lennard-Jones"));
    const std::size_t contacts = toml::get<std::size_t>(toml_value_at(
        global, "ignored_contacts", "[forcefield.global] for Lennard-Jones"));
    MJOLNIR_LOG_INFO("particles connected less than ", bonds,    " bonds are ignored");
    MJOLNIR_LOG_INFO("particles connected less than ", contacts, " contacts are ignored");

    const auto& ps = toml_value_at(global, "parameters",
        "[forcefield.global] for Lennard-Jones potential"
        ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    MJOLNIR_LOG_INFO("{{{");
    for(const auto& param : ps)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.global.parameters]] for Lennard-Jones"));
        if(params.size() <= idx)
        {
            const std::pair<real_type, real_type> dummy{0., 0.};
            params.resize(idx+1, dummy);
        }

        const auto sigma   = toml::get<real_type>(toml_value_at(tab, "sigma",
            "element of [[forcefield.global.parameters]] for Lennard-Jones"));
        const auto epsilon = toml::get<real_type>(toml_value_at(tab, "epsilon",
            "element of [[forcefield.global.parameters]] for Lennard-Jones"));
        MJOLNIR_LOG_INFO("idx = ", idx);
        MJOLNIR_LOG_INFO("sgm = ", sigma);
        MJOLNIR_LOG_INFO("eps = ", epsilon);

        params.at(idx) = std::make_pair(sigma, epsilon);
    }
    MJOLNIR_LOG_INFO("}}}");

    return LennardJonesPotential<traitsT, ignoreT>(
            std::move(params), bonds, contacts);
}

template<typename traitsT, typename ignoreT>
DebyeHuckelPotential<traitsT, ignoreT>
read_debye_huckel_potential(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_debye_huckel_potential(), 0);
    typedef typename traitsT::real_type real_type;

    const std::size_t bonds = toml::get<std::size_t>(toml_value_at(
        global, "ignored_bonds",    "[forcefield.global] for Debye-Huckel"));
    const std::size_t contacts = toml::get<std::size_t>(toml_value_at(
        global, "ignored_contacts", "[forcefield.global] for Debye-Huckel"));
    MJOLNIR_LOG_INFO("particles connected less than ", bonds,    " bonds are ignored");
    MJOLNIR_LOG_INFO("particles connected less than ", contacts, " contacts are ignored");

    const auto& ps = toml_value_at(global, "parameters",
        "[forcefield.global] for Debye-Huckel"
        ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());
    MJOLNIR_LOG_INFO("{{{");
    for(const auto& param : ps)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.global.parameters]] for Debye-Huckel"));
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        const auto charge = toml::get<real_type>(toml_value_at(tab, "charge",
            "element of [[forcefield.global.parameters]] for Debye-Huckel"));
        MJOLNIR_LOG_INFO("idx    = ", idx);
        MJOLNIR_LOG_INFO("charge = ", charge);
        params.at(idx) = charge;
    }
    MJOLNIR_LOG_INFO("}}}");

    return DebyeHuckelPotential<traitsT, ignoreT>(
            std::move(params), bonds, contacts);
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

    const auto thickness = toml::get<real_type>(toml_value_at(external,
        "thickness", "[forcefield.external] for ImplicitMembrane"));
    const auto magnitude = toml::get<real_type>(toml_value_at(external,
        "interaction_magnitude", "[forcefield.external] for ImplicitMembrane"));
    const auto bend = toml::get<real_type>(toml_value_at(external,
        "bend", "[forcefield.external] for ImplicitMembrane"));
    MJOLNIR_LOG_INFO("thickness = ", thickness);
    MJOLNIR_LOG_INFO("magnitude = ", magnitude);
    MJOLNIR_LOG_INFO("bend      = ", bend     );

    const auto& ps = toml_value_at(external, "parameters",
            "[forcefield.external] for ImplicitMembrane"
            ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<real_type> params;
    params.reserve(ps.size());
    MJOLNIR_LOG_INFO("{{{");
    for(const auto& param : ps)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.external.parameters]] for ImplicitMembrane"
            ));

        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        const auto h = toml::get<real_type>(toml_value_at(tab, "hydrophobicity",
            "element of [[forcefield.external.parameters]] for ImplicitMembrane"
            ));
        MJOLNIR_LOG_INFO("idx = ", idx);
        MJOLNIR_LOG_INFO("h   = ", h);
        params.at(idx) = h;
    }
    MJOLNIR_LOG_INFO("}}}");

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

    const auto& ps = toml_value_at(external, "parameters",
            "[forcefield.external] for Lennard-Jones Wall"
            ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    MJOLNIR_LOG_INFO("{{{");
    for(const auto& param : ps)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.external.parameters]] for LennardJonesWall"
            ));
        if(params.size() <= idx)
        {
            params.resize(idx+1, std::make_pair(real_type(0.0), real_type(0.0)));
        }

        const auto s = toml::get<real_type>(toml_value_at(tab, "sigma",
            "element of [[forcefield.external.parameters]] for LennardJonesWall"
            ));
        const auto e = toml::get<real_type>(toml_value_at(tab, "epsilon",
            "element of [[forcefield.external.parameters]] for LennardJonesWall"
            ));
        MJOLNIR_LOG_INFO("idx     = ", idx);
        MJOLNIR_LOG_INFO("sigma   = ", s);
        MJOLNIR_LOG_INFO("epsilon = ", e);
        params.at(idx) = std::make_pair(s, e);
    }
    MJOLNIR_LOG_INFO("}}}");
    return LennardJonesWallPotential<traitsT>(std::move(params));
}

template<typename traitsT>
ExcludedVolumeWallPotential<traitsT>
read_excluded_volume_wall_potential(const toml::Table& external)
{
    typedef typename traitsT::real_type real_type;
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_excluded_volume_wall_potential(), 0);

    const auto& ps = toml_value_at(external, "parameters",
            "[forcefield.external] for Excluded Volume Wall"
            ).cast<toml::value_t::Array>();
    MJOLNIR_LOG_INFO("number of parameters = ", ps.size());

    const real_type eps = toml::get<real_type>(toml_value_at(
        external, "epsilon", "[forcefield.external] for ExcludedVolume"));
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        MJOLNIR_SCOPE(for each parameters, 1);
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.external.parameters]] for ExcludedVolumeWall"
            ));
        if(params.size() <= idx)
        {
            params.resize(idx+1, real_type(0.0));
        }

        const auto s = toml::get<real_type>(toml_value_at(tab, "sigma",
            "element of [[forcefield.external.parameters]] for ExcludedVolumeWall"
            ));
        MJOLNIR_LOG_INFO("idx   = ", idx);
        MJOLNIR_LOG_INFO("sigma = ", s);
        params.at(idx) = s;
    }
    return ExcludedVolumeWallPotential<traitsT>(eps, std::move(params));
}


} // mjolnir
#endif // MJOLNIR_READ_POTENTIAL
