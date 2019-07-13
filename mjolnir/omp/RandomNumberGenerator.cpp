#include <mjolnir/omp/RandomNumberGenerator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class RandomNumberGenerator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class RandomNumberGenerator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class RandomNumberGenerator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class RandomNumberGenerator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
