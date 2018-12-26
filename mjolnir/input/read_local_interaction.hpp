#ifndef MJOLNIR_READ_LOCAL_INTERACTION_HPP
#define MJOLNIR_READ_LOCAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/BondLengthInteraction.hpp>
#include <mjolnir/interaction/BondAngleInteraction.hpp>
#include <mjolnir/interaction/DihedralAngleInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_potential.hpp>
#include <memory>

namespace mjolnir
{

// ----------------------------------------------------------------------------
// local interaction
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_length_interaction(const std::string& kind, const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_bond_length_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[forcefield.local]");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "Go1012Contact")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is 10-12 Go contact.");
        using potentialT = Go1012ContactPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as BondLengthInteraction: ", potential);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_angle_interaction(
    const typename LocalInteractionBase<traitsT>::connection_kind_type kind,
    const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_bond_angle_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[[forcefields.local]]");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "FlexibleLocalAngle")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Flexible Local Angle.");
        using potentialT = FlexibleLocalAnglePotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as BondAngleInteraction: " + potential);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_dihedral_angle_interaction(
    const typename LocalInteractionBase<traitsT>::connection_kind_type kind,
    const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_dihedral_angle_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[[forcefields.local]]");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "ClementiDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Clementi-Go's dihedral.");
        using potentialT = ClementiDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian.");
        using potentialT = AngularGaussianPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "FlexibleLocalDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Flexible Local Dihedral.");
        using potentialT = FlexibleLocalDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    // XXX generalization of this feature is too difficult (technically, it's
    // not so difficult, but practically, it makes the code messy...).
    else if(potential == "Gaussian+FlexibleLocalDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian + FlexibleLocalDihedral.");
        using potential_1_T = GaussianPotential<real_type>;
        using potential_2_T = FlexibleLocalDihedralPotential<real_type>;
        using potentialT    =
            SumLocalPotential<real_type, potential_1_T, potential_2_T>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(kind,
            read_local_potentials<4, real_type, potential_1_T, potential_2_T>(
                local, "Gaussian", "FlexibleLocalDihedral"));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as DihedralAngleInteraction: " + potential);
    }
}

} // mjolnir
#endif// MJOLNIR_READ_LOCAL_INTERACTION
