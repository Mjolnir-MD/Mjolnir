#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
    
template class GlobalStoichiometricInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class GlobalStoichiometricInteraction<SimulatorTraits<float , UnlimitedBoundary>       >;
template class GlobalStoichiometricInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GlobalStoichiometricInteraction<SimulatorTraits<float , CuboidalPeriodicBoundary>>;

} // mjolnir
