#include <mjolnir/core/ObserverContainer.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ObserverContainer<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ObserverContainer<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ObserverContainer<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ObserverContainer<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
