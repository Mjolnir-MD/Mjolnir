#include <mjolnir/forcefield/global/GlobalPairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// // D-H
// template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
// template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
// template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
// template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

} // mjolnir
