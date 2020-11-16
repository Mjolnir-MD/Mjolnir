#include <mjolnir/core/gBAOABLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class gBAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class gBAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class gBAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class gBAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
