#include <mjolnir/forcefield/hybrid/DynamicHybridForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class DynamicHybridForceField<SimulatorTraits<double, UnlimitedBoundary       >>;
template class DynamicHybridForceField<SimulatorTraits<float,  UnlimitedBoundary       >>;
template class DynamicHybridForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class DynamicHybridForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
