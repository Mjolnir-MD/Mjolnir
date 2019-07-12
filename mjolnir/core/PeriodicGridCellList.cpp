#include <mjolnir/core/PeriodicGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<SimulatorTraits<double, UnlimitedBoundary>, empty_t>;
template class PeriodicGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, empty_t>;
template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, empty_t>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, empty_t>;

template class PeriodicGridCellList<SimulatorTraits<double, UnlimitedBoundary>, double>;
template class PeriodicGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, float >;
template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, double>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, float >;

template class PeriodicGridCellList<SimulatorTraits<double, UnlimitedBoundary>, std::pair<double, double>>;
template class PeriodicGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, std::pair<float , float >>;
template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, std::pair<double, double>>;
template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, std::pair<float , float >>;
} // mjolnir
