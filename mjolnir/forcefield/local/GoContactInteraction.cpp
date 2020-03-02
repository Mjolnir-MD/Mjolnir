#include <mjolnir/interaction/local/GoContactInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GoContactPotential<double>>;
template class ContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GoContactPotential<float> >;
template class ContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
template class ContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;
} // mjolnir
