#include <mjolnir/core/ExternalForceInteractionBase.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>>;
template class ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>>;
template class ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
