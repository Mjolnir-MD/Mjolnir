#include <mjolnir/potential/global/LennardJonesPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
