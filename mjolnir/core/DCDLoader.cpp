#include <mjolnir/core/DCDLoader.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class DCDLoader<SimulatorTraits<double, UnlimitedBoundary>>;
template class DCDLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
template class DCDLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class DCDLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
