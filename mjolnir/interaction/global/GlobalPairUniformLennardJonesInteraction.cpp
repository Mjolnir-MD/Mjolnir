#include <mjolnir/interaction/global/GlobalPairUniformLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// Uniform L-J
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> >;

} // mjolnir
