#include <mjolnir/core/MsgPackLoader.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MsgPackLoader<SimulatorTraits<double, UnlimitedBoundary>>;
template class MsgPackLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
template class MsgPackLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class MsgPackLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
