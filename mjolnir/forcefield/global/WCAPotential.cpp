#include <mjolnir/forcefield/global/WCAPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class WCAPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class WCAPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class WCAPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class WCAPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
