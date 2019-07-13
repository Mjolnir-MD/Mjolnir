#include <mjolnir/core/EnergyObserver.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class EnergyObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
template class EnergyObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class EnergyObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class EnergyObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
