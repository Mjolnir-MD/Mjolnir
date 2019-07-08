#include <mjolnir/core/BAOABLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
