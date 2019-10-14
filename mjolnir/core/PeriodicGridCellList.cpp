#include <mjolnir/core/PeriodicGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>,         true>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >,         true>;

template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>,      true>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >,      true>;

template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>,        true>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >,        true>;

template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, true>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >, true>;
} // mjolnir
