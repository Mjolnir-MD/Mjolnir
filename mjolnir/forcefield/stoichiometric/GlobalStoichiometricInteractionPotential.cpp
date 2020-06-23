#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteractionPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalStoichiometricInteractionPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class GlobalStoichiometricInteractionPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class GlobalStoichiometricInteractionPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GlobalStoichiometricInteractionPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
