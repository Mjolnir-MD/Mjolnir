#include <mjolnir/core/ExternalForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ExternalForceField<SimulatorTraits<double, UnlimitedBoundary>>;
template class ExternalForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
template class ExternalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ExternalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
