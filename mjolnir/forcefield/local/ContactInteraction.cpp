#include <mjolnir/forcefield/local/ContactInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// go-contact
template class ContactInteraction<SimulatorTraits<double, UnlimitedBoundary>, GoContactPotential<double>>;
template class ContactInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GoContactPotential<float> >;
template class ContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
template class ContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;

// gaussian
template class ContactInteraction<SimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
template class ContactInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
template class ContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
template class ContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

} // mjolnir
