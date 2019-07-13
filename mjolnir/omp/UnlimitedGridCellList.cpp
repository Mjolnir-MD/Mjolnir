#include <mjolnir/omp/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, empty_t>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, empty_t>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, empty_t>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, empty_t>;

template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, double>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, float >;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, double>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, float >;

template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, std::pair<double, double>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, std::pair<float , float >>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, std::pair<double, double>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, std::pair<float , float >>;
} // mjolnir
