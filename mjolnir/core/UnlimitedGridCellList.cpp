#include <mjolnir/core/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;
} // mjolnir
