#include <mjolnir/forcefield/hybrid/HybridForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class HybridForceField<SimulatorTraits<double, UnlimitedBoundary       >>;
template class HybridForceField<SimulatorTraits<float,  UnlimitedBoundary       >>;
template class HybridForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class HybridForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
