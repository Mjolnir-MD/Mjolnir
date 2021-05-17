#include <mjolnir/forcefield/global/TabulatedWCAPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class TabulatedWCAPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class TabulatedWCAPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class TabulatedWCAPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class TabulatedWCAPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
