#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class GlobalStoichiometricInteraction<SimulatorTraits<double, UnlimitedBoundary>,GlobalStoichiometricInteractionPotential<double>>;
template class GlobalStoichiometricInteraction<SimulatorTraits<float , UnlimitedBoundary>,GlobalStoichiometricInteractionPotential<float >>;
template class GlobalStoichiometricInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GlobalStoichiometricInteractionPotential<double>>;
template class GlobalStoichiometricInteraction<SimulatorTraits<float , CuboidalPeriodicBoundary>, GlobalStoichiometricInteractionPotential<float >>;

} // mjolnir
