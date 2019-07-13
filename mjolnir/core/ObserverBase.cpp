#include <mjolnir/core/ObserverBase.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ObserverBase<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ObserverBase<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ObserverBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ObserverBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
