#include <mjolnir/core/GlobalForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalForceField<SimulatorTraits<double, UnlimitedBoundary>>;
template class GlobalForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
template class GlobalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GlobalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
