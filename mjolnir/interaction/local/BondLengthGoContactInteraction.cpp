#include <mjolnir/interaction/local/BondLengthGoContactInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class BondLengthInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GoContactPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GoContactPotential<float> >;
template class BondLengthInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;

} // mjolnir
