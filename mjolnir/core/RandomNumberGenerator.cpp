#include <mjolnir/core/RandomNumberGenerator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class RandomNumberGenerator<SimulatorTraits<double, UnlimitedBoundary>>;
template class RandomNumberGenerator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class RandomNumberGenerator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class RandomNumberGenerator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
