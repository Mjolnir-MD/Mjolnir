#include <mjolnir/omp/UnderdampedLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
