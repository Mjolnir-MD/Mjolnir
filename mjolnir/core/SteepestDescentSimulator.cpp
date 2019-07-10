#include <mjolnir/core/SteepestDescentSimulator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// BAOAB
template class SteepestDescentSimulator<SimulatorTraits<double, UnlimitedBoundary>       >;
template class SteepestDescentSimulator<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class SteepestDescentSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class SteepestDescentSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
