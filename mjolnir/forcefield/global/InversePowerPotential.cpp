#include <mjolnir/forcefield/global/InversePowerPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class InversePowerParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class InversePowerParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class InversePowerParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class InversePowerParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
}// mjolnir
