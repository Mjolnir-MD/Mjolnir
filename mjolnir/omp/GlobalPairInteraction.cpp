#include <mjolnir/omp/GlobalPairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// D-H
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
