#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class GlobalStoichiometricInteraction<SimulatorTraits<double, UnlimitedBoundary>,StoichiometricUniformCubicPanPotential<double>>;
template class GlobalStoichiometricInteraction<SimulatorTraits<float , UnlimitedBoundary>,StoichiometricUniformCubicPanPotential<float >>;
template class GlobalStoichiometricInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, StoichiometricUniformCubicPanPotential<double>>;
template class GlobalStoichiometricInteraction<SimulatorTraits<float , CuboidalPeriodicBoundary>, StoichiometricUniformCubicPanPotential<float >>;

} // mjolnir
