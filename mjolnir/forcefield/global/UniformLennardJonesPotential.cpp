#include <mjolnir/forcefield/global/UniformLennardJonesPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
