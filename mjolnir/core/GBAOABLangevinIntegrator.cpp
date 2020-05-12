#include <mjolnir/core/GBAOABLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GBAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class GBAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class GBAOABLangevinIntegrator<SimulatorTraits<double, CuboidalBoundary>>;
template class GBAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalBoundary>>;
} // mjolnir
