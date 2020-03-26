#include <mjolnir/forcefield/global/GlobalPairUniformLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// Uniform L-J
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

} // mjolnir
