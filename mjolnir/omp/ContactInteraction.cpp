#include <mjolnir/omp/ContactInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ContactInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        GaussianPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        GaussianPotential<float> >;
template class ContactInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;
} // mjolnir
