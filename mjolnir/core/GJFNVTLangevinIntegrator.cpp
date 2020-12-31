#include <mjolnir/core/GJFNVTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GJFNVTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class GJFNVTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class GJFNVTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GJFNVTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
