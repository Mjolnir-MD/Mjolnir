#include <mjolnir/omp/GlobalStoichiometricInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, StoichiometricUniformCubicPanPotential<double>>;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float , UnlimitedBoundary>, StoichiometricUniformCubicPanPotential<float >>;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, StoichiometricUniformCubicPanPotential<double>>;
template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float , CuboidalPeriodicBoundary>, StoichiometricUniformCubicPanPotential<float >>;

} // mjolnir
