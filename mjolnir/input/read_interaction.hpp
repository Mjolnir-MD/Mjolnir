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
#include <mjolnir/input/get_toml_value.hpp>
#include <mjolnir/input/read_potential.hpp>
#include <mjolnir/input/read_spatial_partition.hpp>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_length_interaction(const toml::Table& local)
{
    const auto potential = toml::get<std::string>(
            toml_value_at(local, "potential", "[forcefield.local]"));
    if(potential == "Harmonic")
    {
        return make_unique<BondLengthInteraction<
			       traitsT, HarmonicPotential<traitsT>>>(
				   read_harmonic_potential<traitsT, 2>(local));
    }
    else if(potential == "Go1012Contact")
    {
        return make_unique<BondLengthInteraction<
			       traitsT, Go1012ContactPotential<traitsT>>>(
				   read_go1012_contact_potential<traitsT, 2>(local));
    }
    else if(potential == "AICG2PlusAngle")
    {
        return make_unique<BondLengthInteraction<
			       traitsT, GaussianPotential<traitsT>>>(
				   read_gaussian_potential<traitsT, 2>(local));
    }
    else
    {
        throw std::runtime_error("invalid length potential: " + potential);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_angle_interaction(const toml::Table& local)
{
    const auto potential = toml::get<std::string>(
            toml_value_at(local, "potential", "[[forcefield.local]]"));
    if(potential == "Harmonic")
    {
        return make_unique<BondAngleInteraction<
			       traitsT, HarmonicPotential<traitsT>>>(
				   read_harmonic_potential<traitsT, 3>(local));
    }
    else if(potential == "FlexibleLocalAngle")
    {
        return make_unique<BondAngleInteraction<
			       traitsT, FlexibleLocalAnglePotential<traitsT>>>(
				   read_flexible_local_angle_potential<traitsT, 3>(local));
    }
    else
    {
        throw std::runtime_error("invalid angle potential: " + potential);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_dihedral_angle_interaction(const toml::Table& local)
{
    const auto potential = toml::get<std::string>(
            toml_value_at(local, "potential", "[forcefield.local]"));
    if(potential == "Harmonic")
    {
        return make_unique<DihedralAngleInteraction<
			       traitsT, HarmonicPotential<traitsT>>>(
				   read_harmonic_potential<traitsT, 4>(local));
    }
    else if(potential == "ClementiDihedral")
    {
        return make_unique<DihedralAngleInteraction<
			       traitsT, ClementiDihedralPotential<traitsT>>>(
				   read_clementi_dihedral_potential<traitsT, 4>(local));
    }
    else if(potential == "AICG2PlusDihedral")
    {
        return make_unique<DihedralAngleInteraction<
			       traitsT, GaussianPotential<traitsT>>>(
				   read_gaussian_potential<traitsT, 4>(local));
    }
    else if(potential == "FlexibleLocalDihedral")
    {
        return make_unique<DihedralAngleInteraction<
			       traitsT, FlexibleLocalDihedralPotential<traitsT>>>(
				   read_flexible_local_dihedral_potential<traitsT, 4>(local));
    }
    else
    {
        throw std::runtime_error("invalid dihedral potential: " + potential);
    }
}

template<typename traitsT, template<typename> class ignoreT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_distance_interaction(const toml::Table& global)
{
    const auto potential = toml::get<std::string>(
            toml_value_at(global, "potential", "[forcefield.local]"));
    if(potential == "ExcludedVolume")
    {
        return read_spatial_partition_for_distance<
            traitsT, ExcludedVolumePotential<traitsT, ignoreT>>(global,
                read_excluded_volume_potential<traitsT, ignoreT>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        return read_spatial_partition_for_distance<
            traitsT, DebyeHuckelPotential<traitsT, ignoreT>>(global,
		        read_debye_huckel_potential<traitsT, ignoreT>(global));
    }
    else if(potential == "LennardJones")
    {
        return read_spatial_partition_for_distance<
            traitsT, LennardJonesPotential<traitsT, ignoreT>>(global,
                read_lennard_jones_potential<traitsT, ignoreT>(global));
    }
    else
    {
        throw std::runtime_error("invalid distance potential: " + potential);
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

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction(const toml::Table& local)
{
    const auto interaction = toml::get<std::string>(
            toml_value_at(local, "interaction", "[forcefields.local]"));
    if(interaction == "BondLength")
    {
        return read_bond_length_interaction<traitsT>(local);
    }
    else if(interaction == "BondAngle")
    {
        return read_bond_angle_interaction<traitsT>(local);
    }
    else if(interaction == "DihedralAngle")
    {
        return read_dihedral_angle_interaction<traitsT>(local);
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
    const auto ignored_group = toml::get<std::string>(
            toml_value_at(global, "ignored_group", "[forcefields.global]"));
    if(interaction == "Distance")
    {
        if(ignored_group == "Nothing")
        {
            return read_global_distance_interaction<
                traitsT, IgnoreNothing>(global);
        }
        else if(ignored_group == "Self")
        {
            return read_global_distance_interaction<
                traitsT, IgnoreSelf>(global);
        }
        else if(ignored_group == "Others")
        {
            return read_global_distance_interaction<
                traitsT, IgnoreOthers>(global);
        }
        else
        {
            throw std::runtime_error(
                    "invalid `ignored_group`: " + ignored_group);
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
