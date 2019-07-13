#include <mjolnir/omp/BondAngleInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// harmonic
template class BondAngleInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, HarmonicPotential<double>>;
template class BondAngleInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, HarmonicPotential<float> >;
template class BondAngleInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
template class BondAngleInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// gaussian
template class BondAngleInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
template class BondAngleInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
template class BondAngleInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
template class BondAngleInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

// FLP angle
template class BondAngleInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, FlexibleLocalAnglePotential<double>>;
template class BondAngleInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, FlexibleLocalAnglePotential<float> >;
template class BondAngleInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<double>>;
template class BondAngleInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<float> >;

} // mjolnir
