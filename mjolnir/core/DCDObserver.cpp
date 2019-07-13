#include <mjolnir/core/DCDObserver.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class DCDObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
template class DCDObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class DCDObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class DCDObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
