#include <mjolnir/omp/GlobalPairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// D-H
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float >>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >>;
} // mjolnir
