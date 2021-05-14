#include <mjolnir/forcefield/global/UniformCubicFunctionPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error  "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UniformCubicFunctionPotential<SimulatorTraits<double, UnlimitedBoundary>>;
template class UniformCubicFunctionPotential<SimulatorTraits<float,  UnlimitedBoundary>>;
template class UniformCubicFunctionPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class UniformCubicFunctionPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
