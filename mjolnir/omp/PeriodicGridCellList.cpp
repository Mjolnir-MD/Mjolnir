#include <mjolnir/omp/PeriodicGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >>;
} // mjolnir
