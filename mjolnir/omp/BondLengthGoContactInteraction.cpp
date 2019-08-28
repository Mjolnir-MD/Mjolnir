#include <mjolnir/omp/BondLengthGoContactInteraction.hpp>

namespace mjolnir
{
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary       >, GoContactPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary       >, GoContactPotential<float> >;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;
} // mjolnir
