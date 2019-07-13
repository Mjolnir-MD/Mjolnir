#include <mjolnir/core/XYZObserver.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class XYZObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
template class XYZObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class XYZObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class XYZObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
