#include <mjolnir/omp/GlobalPairUniformLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> >;
} // mjolnir
