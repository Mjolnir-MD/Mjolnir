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

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<float >>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, InversePowerPotential<double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, InversePowerPotential<float >>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, WCAPotential<double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, WCAPotential<float >>;
}

