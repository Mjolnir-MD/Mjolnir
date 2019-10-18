#include <mjolnir/core/PeriodicGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >>;

template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >>;

template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >>;
} // mjolnir
