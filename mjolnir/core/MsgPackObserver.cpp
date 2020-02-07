#include <mjolnir/core/MsgPackObserver.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MsgPackObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
template class MsgPackObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class MsgPackObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class MsgPackObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
