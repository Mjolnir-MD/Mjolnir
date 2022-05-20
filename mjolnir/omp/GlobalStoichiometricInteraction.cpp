#include <mjolnir/omp/GlobalStoichiometricInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, GlobalStoichiometricInteractionPotential<double>>;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float , UnlimitedBoundary>, GlobalStoichiometricInteractionPotential<float >>;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GlobalStoichiometricInteractionPotential<double>>;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float , CuboidalPeriodicBoundary>, GlobalStoichiometricInteractionPotential<float >>;

} // mjolnir
