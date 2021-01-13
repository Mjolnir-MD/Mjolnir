#include <mjolnir/omp/GlobalStoichiometricInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
    
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float , UnlimitedBoundary>       >;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float , CuboidalPeriodicBoundary>>;

} // mjolnir
