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

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, InversePowerPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, InversePowerPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, WCAPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, WCAPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
