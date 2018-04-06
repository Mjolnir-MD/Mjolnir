#ifndef MJOLNIR_READ_INTERACTION
#define MJOLNIR_READ_INTERACTION
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/ZaxisExternalForceInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/input/get_toml_value.hpp>
#include <mjolnir/input/read_potential.hpp>
#include <mjolnir/input/read_spatial_partition.hpp>
#include <memory>

namespace mjolnir
{

// ----------------------------------------------------------------------------
// local interaction
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_length_interaction(
    const typename LocalInteractionBase<traitsT>::connection_kind_type kind,
    const toml::Table& local)
{
    const auto potential = toml::get<std::string>(
            toml_value_at(local, "potential", "[forcefield.local]"));

    if(potential == "Harmonic")
    {
        using potential_t = HarmonicPotential<traitsT>;

        return make_unique<BondLengthInteraction<traitsT, potential_t>>(
            kind, read_harmonic_potential<traitsT, 2>(local));
    }
    else if(potential == "Go1012Contact")
    {
        using potential_t = Go1012ContactPotential<traitsT>;

        return make_unique<BondLengthInteraction<traitsT, potential_t>>(
            kind, read_go1012_contact_potential<traitsT, 2>(local));
    }
    else if(potential == "AICG2PlusAngle")
    {
        using potential_t = GaussianPotential<traitsT>;

        return make_unique<BondLengthInteraction<traitsT, potential_t>>(
            kind, read_gaussian_potential<traitsT, 2>(local));
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
    const auto potential = toml::get<std::string>(
            toml_value_at(local, "potential", "[[forcefield.local]]"));
    if(potential == "Harmonic")
    {
        using potential_t = HarmonicPotential<traitsT>;

        return make_unique<BondAngleInteraction<traitsT, potential_t>>(
            kind, read_harmonic_potential<traitsT, 3>(local));
    }
    else if(potential == "FlexibleLocalAngle")
    {
        using potential_t = FlexibleLocalAnglePotential<traitsT>;

        return make_unique<BondAngleInteraction<traitsT, potential_t>>(
            kind, read_flexible_local_angle_potential<traitsT, 3>(local));
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
    const auto potential = toml::get<std::string>(
            toml_value_at(local, "potential", "[forcefield.local]"));
    if(potential == "Harmonic")
    {
        using potential_t = HarmonicPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potential_t>>(
            kind, read_harmonic_potential<traitsT, 4>(local));
    }
    else if(potential == "ClementiDihedral")
    {
        using potential_t = ClementiDihedralPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potential_t>>(
            kind, read_clementi_dihedral_potential<traitsT, 4>(local));
    }
    else if(potential == "AICG2PlusDihedral")
    {
        using potential_t = GaussianPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potential_t>>(
            kind, read_gaussian_potential<traitsT, 4>(local));
    }
    else if(potential == "FlexibleLocalDihedral")
    {
        using potential_t = FlexibleLocalDihedralPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potential_t>>(
            kind, read_flexible_local_dihedral_potential<traitsT, 4>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as DihedralAngleInteraction: " + potential);
    }
}

// ----------------------------------------------------------------------------
// global interaction
// ----------------------------------------------------------------------------

template<typename traitsT, typename ignoreT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_distance_interaction(const toml::Table& global)
{
    const auto potential = toml::get<std::string>(
            toml_value_at(global, "potential", "[forcefield.local]"));
    if(potential == "ExcludedVolume")
    {
        using potential_t = ExcludedVolumePotential<traitsT, ignoreT>;

        return read_spatial_partition_for_distance<traitsT, potential_t>(
            global, read_excluded_volume_potential<traitsT, ignoreT>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        using potential_t = DebyeHuckelPotential<traitsT, ignoreT>;

        return read_spatial_partition_for_distance<traitsT, potential_t>(
            global, read_debye_huckel_potential<traitsT, ignoreT>(global));
    }
    else if(potential == "LennardJones")
    {
        using potential_t = LennardJonesPotential<traitsT, ignoreT>;

        return read_spatial_partition_for_distance<traitsT, potential_t>(
            global, read_lennard_jones_potential<traitsT, ignoreT>(global));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as GlobalDistanceInteraction: " + potential);
    }
}

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_zaxis_external_force_interaction(const toml::Table& global)
{
    const auto potential = toml::get<std::string>(
            toml_value_at(global, "potential", "[forcefield.local]"));
    if(potential == "ImplicitMembrane")
    {
        return read_spatial_partition_for_implicit_membrane<
            traitsT, ImplicitMembranePotential<traitsT>>(
                global, read_implicit_membrane_potential<traitsT>(global));
    }
    else
    {
        throw std::runtime_error("invalid distance potential: " + potential);
    }
}

// ----------------------------------------------------------------------------
// general read_(local|global)_interaction function
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction(const toml::Table& local)
{
    const auto interaction = toml::get<std::string>(
            toml_value_at(local, "interaction", "[forcefields.local]"));

    // topology stuff
    using connection_kind_type =
        typename LocalInteractionBase<traitsT>::connection_kind_type;
    const auto connection = toml::get<std::string>(
            toml_value_at(local, "topology", "[forcefield.local]"));

    connection_kind_type kind;
    if     (connection == "bond")    {kind = connection_kind_type::bond;}
    else if(connection == "contact") {kind = connection_kind_type::contact;}
    else if(connection == "none")    {kind = connection_kind_type::none;}
    else {throw std::runtime_error("invalid connection type: " + connection);}

    if(interaction == "BondLength")
    {
        return read_bond_length_interaction<traitsT>(kind, local);
    }
    else if(interaction == "BondAngle")
    {
        return read_bond_angle_interaction<traitsT>(kind, local);
    }
    else if(interaction == "DihedralAngle")
    {
        return read_dihedral_angle_interaction<traitsT>(kind, local);
    }
    else
    {
        throw std::runtime_error(
                "invalid local interaction type: " + interaction);
    }
}

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_interaction(const toml::Table& global)
{
    const auto interaction = toml::get<std::string>(
            toml_value_at(global, "interaction", "[forcefields.global]"));
    const auto ignored_chain = toml::get<std::string>(
            toml_value_at(global, "ignored_chain", "[forcefields.global]"));

    if(interaction == "Distance")
    {
        if(ignored_chain == "Nothing")
        {
            return read_global_distance_interaction<traitsT, IgnoreNothing>(global);
        }
        else if(ignored_chain == "Self")
        {
            return read_global_distance_interaction<traitsT, IgnoreSelf>(global);
        }
        else if(ignored_chain == "Others")
        {
            return read_global_distance_interaction<traitsT, IgnoreOthers>(global);
        }
        else
        {
            throw std::runtime_error("invalid `ignored_chain`: " + ignored_chain);
        }
    }
    else if(interaction == "External")
    {
        return read_zaxis_external_force_interaction<traitsT>(global);
    }
    else
    {
        throw std::runtime_error(
                "invalid global interaction type: " + interaction);
    }
}

} // mjolnir
#endif// MJOLNIR_READ_INTERACTION
