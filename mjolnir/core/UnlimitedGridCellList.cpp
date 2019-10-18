#include <mjolnir/core/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<double>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<float >>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<double>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<float >>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<double>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<float >>;
} // mjolnir
