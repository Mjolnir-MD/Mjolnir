#include <mjolnir/omp/ContactInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// go-contact
template class ContactInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, GoContactPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, GoContactPotential<float> >;
template class ContactInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;

// gaussian
template class ContactInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
template class ContactInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
template class ContactInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

} // mjolnir
