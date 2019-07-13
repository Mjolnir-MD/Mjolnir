#include <mjolnir/omp/PeriodicGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, empty_t>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, empty_t>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, empty_t>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, empty_t>;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, double>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, float >;
template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, double>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, float >;

template class PeriodicGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, std::pair<double, double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, std::pair<float , float >>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, std::pair<double, double>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, std::pair<float , float >>;
} // mjolnir
