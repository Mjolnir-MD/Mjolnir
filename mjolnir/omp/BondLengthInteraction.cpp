#include <mjolnir/omp/BondLengthInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// harmonic
template class BondLengthInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, HarmonicPotential<double>>;
template class BondLengthInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, HarmonicPotential<float> >;
template class BondLengthInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
template class BondLengthInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// gaussian
template class BondLengthInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
template class BondLengthInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
template class BondLengthInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
template class BondLengthInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

} // mjolnir
