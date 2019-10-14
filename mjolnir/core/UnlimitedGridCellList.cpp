#include <mjolnir/core/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<double>,         true>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<float >,         true>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<double>,      true>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<float >,      true>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<double>,        true>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<float >,        true>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<double>, true>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<float >, true>;
} // mjolnir
