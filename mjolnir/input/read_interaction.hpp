#ifndef MJOLNIR_READ_INTERACTION
#define MJOLNIR_READ_INTERACTION
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/GlobalExternalInteractionMU.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <memory>
#include "read_potential.hpp"
#include "read_spatial_partition.hpp"

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_length_interaction(const toml::Table& local)
{
    const auto potential = toml::get<std::string>(local.at("potential"));
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
    const auto potential = toml::get<std::string>(local.at("potential"));
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
    const auto potential = toml::get<std::string>(local.at("potential"));
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

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_distance_interaction(const toml::Table& global)
{
    const auto potential = toml::get<std::string>(global.at("potential"));
    if(potential == "ExcludedVolume")
    {
        return read_spatial_partition_for_distance<
            traitsT, ExcludedVolumePotential<traitsT>>(
                       global, read_excluded_volume_potential<traitsT>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        return read_spatial_partition_for_distance<
            traitsT, DebyeHuckelPotential<traitsT>>(
                       global, read_debye_huckel_potential<traitsT>(global));
    }
    else if(potential == "LennardJones")
    {
        return read_spatial_partition_for_distance<
            traitsT, LennardJonesPotential<traitsT>>(
                       global, read_lennard_jones_potential<traitsT>(global));
    }
    else
    {
        throw std::runtime_error("invalid distance potential: " + potential);
    }
}

  template<typename traitsT>
  std::unique_ptr<GlobalInteractionBase<traitsT>>
  read_global_external_interaction(const toml::Table& global)
  {
    const auto potential = toml::get<std::string>(global.at("potential"));
    if(potential == "ImplicitMembrane")
      {
	//This is for MU.
	return read_spatial_partition_for_externalMU<
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
    const auto interaction = toml::get<std::string>(local.at("interaction"));
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
    const auto interaction = toml::get<std::string>(global.at("interaction"));
    if(interaction == "Distance")
    {
        return read_global_distance_interaction<traitsT>(global);
    }
    else if(interaction == "External")
      {
	return read_global_external_interaction<traitsT>(global);
      }
    else
    {
        throw std::runtime_error(
                "invalid global interaction type: " + interaction);
    }
}

} // mjolnir
#endif// MJOLNIR_READ_INTERACTION
