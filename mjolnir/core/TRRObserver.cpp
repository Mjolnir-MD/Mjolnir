#include <mjolnir/core/TRRObserver.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class TRRObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
template class TRRObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class TRRObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class TRRObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
