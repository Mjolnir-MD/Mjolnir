#include <mjolnir/core/MsgPackSaver.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MsgPackSaver<SimulatorTraits<double, UnlimitedBoundary>       >;
template class MsgPackSaver<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class MsgPackSaver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class MsgPackSaver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
