#include <mjolnir/interaction/local/BondLengthInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// harmonic
template class BondLengthInteraction<SimulatorTraits<double, UnlimitedBoundary>, HarmonicPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  UnlimitedBoundary>, HarmonicPotential<float> >;
template class BondLengthInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// go-contact
template class BondLengthInteraction<SimulatorTraits<double, UnlimitedBoundary>, GoContactPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GoContactPotential<float> >;
template class BondLengthInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;

// gaussian
template class BondLengthInteraction<SimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
template class BondLengthInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
template class BondLengthInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

} // mjolnir
