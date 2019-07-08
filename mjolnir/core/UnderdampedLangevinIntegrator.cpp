#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
