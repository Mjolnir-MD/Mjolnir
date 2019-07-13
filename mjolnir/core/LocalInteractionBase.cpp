#include <mjolnir/core/LocalInteractionBase.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>>;
template class LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>>;
template class LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
