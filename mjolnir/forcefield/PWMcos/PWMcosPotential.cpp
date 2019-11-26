#include <mjolnir/forcefield/PWMcos/PWMcosPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PWMcosPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class PWMcosPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class PWMcosPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class PWMcosPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
