#include <mjolnir/core/GlobalInteractionBase.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>>;
template class GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>>;
template class GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
