#include <mjolnir/omp/PeriodicGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
