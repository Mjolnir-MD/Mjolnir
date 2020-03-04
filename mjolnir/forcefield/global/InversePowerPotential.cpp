#include <mjolnir/forcefield/global/InversePowerPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class InversePowerPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class InversePowerPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class InversePowerPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class InversePowerPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
}// mjolnir
