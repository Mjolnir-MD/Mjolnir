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

template<typename traitsT, typename boundaryT, typename potentialT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction_boundary_potential(const std::string& interaction,
        const toml::Array<toml::Table>& parameters)
{
    if(interaction == "BondLength")
    {
        return make_unique<
            BondLengthInteraction<traitsT, potentialT, boundaryT>>(
                read_local_potential<traitsT, potentialT, 2>(parameters));
    }
    else if(interaction == "BondAngle")
    {
        return make_unique<
            BondAngleInteraction<traitsT, potentialT, boundaryT>>(
                read_local_potential<traitsT, potentialT, 3>(parameters));
    }
    else if(interaction == "DihedralAngle")
    {
        return make_unique<
            DihedralAngleInteraction<traitsT, potentialT, boundaryT>>(
                read_local_potential<traitsT, potentialT, 4>(parameters));
    }
    else
        throw std::runtime_error("unknown potential: " + interaction);
}

template<typename traitsT, typename boundaryT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction_boundary(const std::string& interaction,
    const std::string& potential, const toml::Array<toml::Table>& parameters)
{
    if(potential == "Harmonic")
    {
        return read_local_interaction_boundary_potential<traitsT, boundaryT,
            HarmonicPotential<traitsT>>(interaction, parameters);
    }
    else if(potential == "ClementiDihedral")
    {
        return read_local_interaction_boundary_potential<traitsT, boundaryT,
            ClementiDihedralPotential<traitsT>>(interaction, parameters);
    }
    else if(potential == "Go1012Contact")
    {
        return read_local_interaction_boundary_potential<traitsT, boundaryT,
            Go1012ContactPotential<traitsT>>(interaction, parameters);
    }
    else
        throw std::runtime_error("unknown potential: " + potential);
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction(const std::string& interaction,
        const std::string& potential, const std::string& boundary,
        const toml::Array<toml::Table>& parameters)
{
    if(boundary == "Unlimited")
    {
        return read_local_interaction_boundary<traitsT,
            UnlimitedBoundary<traitsT>>(interaction, potential, parameters);
    }
    else if(boundary == "Periodic")
    {
        return read_local_interaction_boundary<traitsT,
            PeriodicBoundaryXYZ<traitsT>>(interaction, potential, parameters);
    }
    else
        throw std::runtime_error("unknown boundary condition: " + boundary);
}

} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_LOCAL_INTERACTION */
