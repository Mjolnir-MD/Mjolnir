#ifndef MJOLNIR_IO_TOML_READ_LOCAL_INTERACTION
#define MJOLNIR_IO_TOML_READ_LOCAL_INTERACTION
#include "read_local_potential.hpp"
#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 2>>
read_unlimited_bond_length_interaction(
        const std::string& name, const toml::Table& params)
{
    if(name == "Harmonic")
        return make_unique<BondLengthInteraction<traitsT,
                   HarmonicPotential<traitsT>, UnlimitedBoundary<traitsT>>>(
                           read_harmonic_potential<traitsT>(params));
    else if(name == "Go1012Contact")
        return make_unique<BondLengthInteraction<traitsT,
                   Go1012ContactPotential<traitsT>, UnlimitedBoundary<traitsT>>>(
                           read_go_1012_contact_potential<traitsT>(params));
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 2>>
read_periodic_bond_length_interaction(
        const std::string& name, const toml::Table& params)
{
    if(name == "Harmonic")
        return make_unique<BondLengthInteraction<traitsT,
                   HarmonicPotential<traitsT>, PeriodicBoundaryXYZ<traitsT>>>(
                           read_harmonic_potential<traitsT>(params));
    else if(name == "Go1012Contact")
        return make_unique<BondLengthInteraction<traitsT,
                   Go1012ContactPotential<traitsT>, PeriodicBoundaryXYZ<traitsT>>>(
                           read_go_1012_contact_potential<traitsT>(params));
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 3>>
read_unlimited_bond_angle_interaction(
        const std::string& name, const toml::Table& params)
{
    if(name == "Harmonic")
        return make_unique<BondAngleInteraction<traitsT,
                   HarmonicPotential<traitsT>, UnlimitedBoundary<traitsT>>>(
                           read_harmonic_potential<traitsT>(params));
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 3>>
read_periodic_bond_angle_interaction(
        const std::string& name, const toml::Table& params)
{
    if(name == "Harmonic")
        return make_unique<BondAngleInteraction<traitsT,
                   HarmonicPotential<traitsT>, PeriodicBoundaryXYZ<traitsT>>>(
                           read_harmonic_potential<traitsT>(params));
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 4>>
read_unlimited_dihedral_angle_interaction(
        const std::string& name, const toml::Table& params)
{
    if(name == "Harmonic")
        return make_unique<DihedralAngleInteraction<traitsT,
                   HarmonicPotential<traitsT>, UnlimitedBoundary<traitsT>>>(
                           read_harmonic_potential<traitsT>(params));
    else if(name == "ClementiDihedral")
        return make_unique<DihedralAngleInteraction<traitsT,
                   ClementiDihedralPotential<traitsT>, UnlimitedBoundary<traitsT>>>(
                           read_clementi_dihedral_potential<traitsT>(params));
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 4>>
read_periodic_dihedral_angle_interaction(
        const std::string& name, const toml::Table& params)
{
    if(name == "Harmonic")
        return make_unique<DihedralAngleInteraction<traitsT,
                   HarmonicPotential<traitsT>, PeriodicBoundaryXYZ<traitsT>>>(
                           read_harmonic_potential<traitsT>(params));
    else if(name == "ClementiDihedral")
        return make_unique<DihedralAngleInteraction<traitsT,
                   ClementiDihedralPotential<traitsT>, PeriodicBoundaryXYZ<traitsT>>>(
                           read_clementi_dihedral_potential<traitsT>(params));
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 2>>
read_2body_interaction(const std::string& name, const std::string& bound,
        const std::string& potential, const toml::Table& params)
{
    if(name == "BondLength")
    {
        if(bound == "Unlimited")
            return read_unlimited_bond_length_interaction<traitsT>(potential, params);
        else if(bound == "Periodic")
            return read_periodic_bond_length_interaction<traitsT>(potential, params);
        else
            throw std::runtime_error("unknown boundary: " + bound);
    }
    else
    {
        throw std::runtime_error(std::string("unknown interaction: ") + name);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 3>>
read_3body_interaction(const std::string& name, const std::string& bound,
        const std::string& potential, const toml::Table& params)
{
    if(name == "BondAngle")
    {
        if(bound == "Unlimited")
            return read_unlimited_bond_angle_interaction<traitsT>(potential, params);
        else if(bound == "Periodic")
            return read_periodic_bond_angle_interaction<traitsT>(potential, params);
        else
            throw std::runtime_error("unknown boundary: " + bound);
    }
    else
    {
        throw std::runtime_error(std::string("unknown interaction: ") + name);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 4>>
read_4body_interaction(const std::string& name, const std::string& bound,
        const std::string& potential, const toml::Table& params)
{
    if(name == "DihedralAngle")
    {
        if(bound == "Unlimited")
            return read_unlimited_dihedral_angle_interaction<traitsT>(potential, params);
        else if(bound == "Periodic")
            return read_periodic_dihedral_angle_interaction<traitsT>(potential, params);
        else
            throw std::runtime_error("unknown boundary: " + bound);
    }
    else
    {
        throw std::runtime_error(std::string("unknown interaction: ") + name);
    }
}

} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_LOCAL_INTERACTION */
