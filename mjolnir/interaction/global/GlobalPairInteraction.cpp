#include <mjolnir/interaction/global/GlobalPairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// D-H
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float> >;

} // mjolnir
