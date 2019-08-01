#include <mjolnir/core/LocalForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LocalForceField<SimulatorTraits<double, UnlimitedBoundary>>;
template class LocalForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
template class LocalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class LocalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
