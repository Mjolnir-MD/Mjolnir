#ifndef MJOLNIR_IO_TOML_READ_GLOBAL_INTERACTION
#define MJOLNIR_IO_TOML_READ_GLOBAL_INTERACTION
#include "read_spatial_partition.hpp"
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename spaceT,
         typename boundaryT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_interaction_boundary_potential(
        const std::string& interaction, const toml::Table& potent)
{
    if(interaction == "Global")
    {
        auto pot = read_global_potential<traitsT, potentialT>(potent);
        const auto cutoff = pot.max_cutoff_length();
        return make_unique<GlobalDistanceInteraction<traitsT,
            potentialT, spaceT, boundaryT>>(std::move(pot),
                read_spatial_partition<traitsT, spaceT, boundaryT>(
                    potent, cutoff));
    }
    else
        throw std::runtime_error("unknown interaction: " + interaction);
}


template<typename traitsT, typename spaceT, typename boundaryT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_interaction_boundary(const std::string& interaction,
        const std::string& potential, const toml::Table& potent)
{
    if(potential == "ExcludedVolume")
    {
        return read_global_interaction_boundary_potential<traitsT,
                ExcludedVolumePotential<traitsT>, spaceT, boundaryT>(
                    interaction, potent);
    }
    else if(potential == "LennardJones")
    {
        return read_global_interaction_boundary_potential<traitsT,
                LennardJonesPotential<traitsT>, spaceT, boundaryT>(
                    interaction, potent);
    }
    else
        throw std::runtime_error("unknown potential: " + potential);
}

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_interaction(const std::string& interaction,
        const std::string& potential, const std::string& boundary,
        const toml::Table& potent)
{
    if(boundary == "Unlimited")
    {
        return read_global_interaction_boundary<traitsT,
                   UnlimitedGridCellList<traitsT>,
                   UnlimitedBoundary<traitsT>>(interaction, potential, potent);
    }
    else if(boundary == "Periodic")
    {
        return read_global_interaction_boundary<traitsT,
                   PeriodicGridCellList<traitsT>,
                   PeriodicBoundaryXYZ<traitsT>>(interaction, potential, potent);
    }
    else
        throw std::runtime_error("unknown interaction: " + interaction);
}



} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_GLOBAL_INTERACTION */
