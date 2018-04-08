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
#include <mjolnir/input/get_toml_value.hpp>

namespace mjolnir
{

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, HarmonicPotential<traitsT>>>
read_harmonic_potential(const toml::Table& local)
{
    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, HarmonicPotential<traitsT>>;

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for Harmonic potential"
        ).cast<toml::value_t::Array>();

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        using real_type = typename traitsT::real_type;

        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in harmonic potential"));
        auto r0 = toml::get<real_type>(toml_value_at(parameter, "native",
                "element of [[parameters]] in harmonic potential"));
        auto k  = toml::get<real_type>(toml_value_at(parameter, "k",
                "element of [[parameters]] in harmonic potential"));

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
    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, Go1012ContactPotential<traitsT>>;

    const auto& params = toml_value_at(local, "parameters",
            "[forcefield.local] for Go-10-12 potential"
            ).cast<toml::value_t::Array>();

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        using real_type = typename traitsT::real_type;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in Go-10-12 potential"));
        auto r0 = toml::get<real_type>(toml_value_at(parameter, "native",
                "element of [[parameters]] in Go-10-12 potential"));
        auto k  = toml::get<real_type>(toml_value_at(parameter, "k",
                "element of [[parameters]] in Go-10-12 potential"));

        retval.emplace_back(indices, Go1012ContactPotential<traitsT>(k, r0));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, GaussianPotential<traitsT>>>
read_gaussian_potential(const toml::Table& local)
{
    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, GaussianPotential<traitsT>>;

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for Gaussian potential"
        ).cast<toml::value_t::Array>();

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        using real_type = typename traitsT::real_type;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in Gaussian potential"));
        auto native  = toml::get<real_type>(toml_value_at(parameter, "native",
                "element of [[parameters]] in Gaussian potential"));
        auto epsilon = toml::get<real_type>(toml_value_at(parameter, "epsilon",
                "element of [[parameters]] in Gaussian potential"));
        auto w       = toml::get<real_type>(toml_value_at(parameter, "w",
                "element of [[parameters]] in Gaussian potential"));

        retval.emplace_back(
                indices, GaussianPotential<traitsT>(epsilon, w, native));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<
    std::pair<std::array<std::size_t, N>, FlexibleLocalAnglePotential<traitsT>>>
read_flexible_local_angle_potential(const toml::Table& local)
{
    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, FlexibleLocalAnglePotential<traitsT>>;

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for FlexibleLocalAngle potential"
        ).cast<toml::value_t::Array>();

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
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
    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, ClementiDihedralPotential<traitsT>>;

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for ClementiDihedral potential"
        ).cast<toml::value_t::Array>();

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        using real_type  = typename traitsT::real_type;
        const auto& parameter = item.cast<toml::value_t::Table>();

        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
                "element of [[parameters]] in ClementiDihedral potential"));
        auto native = toml::get<real_type>(toml_value_at(parameter, "native",
                "element of [[parameters]] in ClementiDihedral potential"));
        auto k1     = toml::get<real_type>(toml_value_at(parameter, "k1",
                "element of [[parameters]] in ClementiDihedral potential"));
        auto k3     = toml::get<real_type>(toml_value_at(parameter, "k3",
                "element of [[parameters]] in ClementiDihedral potential"));

        retval.emplace_back(
                indices, ClementiDihedralPotential<traitsT>(k1, k3, native));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>,
                      FlexibleLocalDihedralPotential<traitsT>>>
read_flexible_local_dihedral_potential(const toml::Table& local)
{
    using indices_t = std::array<std::size_t, N>;
    using indices_potential_pair_t =
        std::pair<indices_t, FlexibleLocalDihedralPotential<traitsT>>;

    const auto& params = toml_value_at(local, "parameters",
        "[forcefield.local] for FlexibleLocalDihedral potential"
        ).cast<toml::value_t::Array>();
    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        using real_type  = typename traitsT::real_type;
        using table_type = std::array<real_type, 7>;

        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<indices_t>(toml_value_at(parameter, "indices",
            "element of [[parameters]] in FlexibleLocalDihedral potential"));
        auto k    = toml::get<real_type>(toml_value_at(parameter, "k",
            "element of [[parameters]] in FlexibleLocalDihedral potential"));
        auto term = toml::get<table_type>(toml_value_at(parameter, "term",
            "element of [[parameters]] in FlexibleLocalDihedral potential"));

        retval.emplace_back(
                indices, FlexibleLocalDihedralPotential<traitsT>(k, term));
    }
    return retval;
}

template<typename traitsT, typename ignoreT>
ExcludedVolumePotential<traitsT, ignoreT>
read_excluded_volume_potential(const toml::Table& global)
{
    typedef typename traitsT::real_type real_type;

    const std::size_t bonds = toml::get<std::size_t>(toml_value_at(
        global, "ignored_bonds", "[forcefield.global] for ExcludedVolume"));
    const std::size_t contacts = toml::get<std::size_t>(toml_value_at(
        global, "ignored_contacts", "[forcefield.global] for ExcludedVolume"));
    const real_type eps = toml::get<real_type>(toml_value_at(
        global, "epsilon", "[forcefield.global] for ExcludedVolume"));

    const auto& ps = toml_value_at(global, "parameters",
            "[forcefield.global] for ExcludedVolume potential"
            ).cast<toml::value_t::Array>();
    std::vector<real_type> params;
    params.reserve(ps.size());

    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.global.parameters]] for Excluded Volume"));
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        params.at(idx) = toml::get<real_type>(toml_value_at(tab, "sigma",
            "element of [[forcefield.global.parameters]] for Excluded Volume"));
    }

    return ExcludedVolumePotential<traitsT, ignoreT>(
            eps, std::move(params), bonds, contacts);
}

template<typename traitsT, typename ignoreT>
LennardJonesPotential<traitsT, ignoreT>
read_lennard_jones_potential(const toml::Table& global)
{
    typedef typename traitsT::real_type real_type;

    const std::size_t bonds = toml::get<std::size_t>(toml_value_at(
        global, "ignored_bonds", "[forcefield.global] for Lennard-Jones"));
    const std::size_t contacts = toml::get<std::size_t>(toml_value_at(
        global, "ignored_contacts", "[forcefield.global] for Lennard-Jones"));

    const auto& ps = toml_value_at(global, "parameters",
        "[forcefield.global] for Lennard-Jones potential"
        ).cast<toml::value_t::Array>();

    std::vector<std::pair<real_type, real_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
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

        params.at(idx) = std::make_pair(sigma, epsilon);
    }

    return LennardJonesPotential<traitsT, ignoreT>(
            std::move(params), bonds, contacts);
}

template<typename traitsT, typename ignoreT>
DebyeHuckelPotential<traitsT, ignoreT>
read_debye_huckel_potential(const toml::Table& global)
{
    typedef typename traitsT::real_type real_type;

    const std::size_t bonds = toml::get<std::size_t>(toml_value_at(
        global, "ignored_bonds",    "[forcefield.global] for Debye-Huckel"));
    const std::size_t contacts = toml::get<std::size_t>(toml_value_at(
        global, "ignored_contacts", "[forcefield.global] for Debye-Huckel"));

    const auto& ps = toml_value_at(global, "parameters",
        "[forcefield.global] for Debye-Huckel"
        ).cast<toml::value_t::Array>();

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(toml_value_at(tab, "index",
            "element of [[forcefield.global.parameters]] for Debye-Huckel"));
        if(params.size() <= idx)
        {
            params.resize(idx+1, 0.);
        }
        const auto charge = toml::get<real_type>(toml_value_at(tab, "charge",
            "element of [[forcefield.global.parameters]] for Debye-Huckel"));
        params.at(idx) = charge;
    }

    return DebyeHuckelPotential<traitsT, ignoreT>(
            std::move(params), bonds, contacts);
}

template<typename traitsT>
ImplicitMembranePotential<traitsT>
read_implicit_membrane_potential(const toml::Table& external)
{
    typedef typename traitsT::real_type real_type;

    const auto thickness = toml::get<real_type>(toml_value_at(external,
        "thickness", "[forcefield.external] for ImplicitMembrane"));
    const auto magnitude = toml::get<real_type>(toml_value_at(external,
        "interaction_magnitude", "[forcefield.external] for ImplicitMembrane"));
    const auto bend = toml::get<real_type>(toml_value_at(external,
        "bend", "[forcefield.external] for ImplicitMembrane"));

    const auto& ps = toml_value_at(external, "parameters",
            "[forcefield.external] for ImplicitMembrane"
            ).cast<toml::value_t::Array>();

    std::vector<real_type> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
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
        params.at(idx) = h;
    }

    return ImplicitMembranePotential<traitsT>(
            thickness, magnitude, bend, std::move(params));
}

} // mjolnir
#endif // MJOLNIR_READ_POTENTIAL
