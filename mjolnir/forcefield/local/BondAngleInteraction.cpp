#include <mjolnir/interaction/local/BondAngleInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// harmonic
template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, HarmonicPotential<double>>;
template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, HarmonicPotential<float> >;
template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// gaussian
template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

// FLP angle
template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, FlexibleLocalAnglePotential<double>>;
template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, FlexibleLocalAnglePotential<float> >;
template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<double>>;
template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<float> >;

} // mjolnir
