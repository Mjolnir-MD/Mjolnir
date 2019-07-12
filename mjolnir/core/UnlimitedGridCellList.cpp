#include <mjolnir/core/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, empty_t>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, empty_t>;
template class UnlimitedGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, empty_t>;
template class UnlimitedGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, empty_t>;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, double>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, float >;
template class UnlimitedGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, double>;
template class UnlimitedGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, float >;

template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, std::pair<double, double>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, std::pair<float , float >>;
template class UnlimitedGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, std::pair<double, double>>;
template class UnlimitedGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, std::pair<float , float >>;
} // mjolnir
