#include <mjolnir/forcefield/PWMcos/PWMcosPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PWMcosParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class PWMcosParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class PWMcosParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class PWMcosParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
