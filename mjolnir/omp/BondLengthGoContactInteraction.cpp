#include <mjolnir/omp/BondLengthGoContactInteraction.hpp>

namespace mjolnir
{
template class BondLengthInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary       >, GoContactPotential<double>>;
template class BondLengthInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary       >, GoContactPotential<float> >;
template class BondLengthInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
template class BondLengthInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;
} // mjolnir
