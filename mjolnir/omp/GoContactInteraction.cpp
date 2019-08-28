#include <mjolnir/omp/GoContactInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ContactInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary       >, GoContactPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary       >, GoContactPotential<float> >;
template class ContactInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;
} // mjolnir
