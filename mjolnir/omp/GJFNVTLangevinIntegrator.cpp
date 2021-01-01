#include <mjolnir/omp/GJFNVTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GJFNVTLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class GJFNVTLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class GJFNVTLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GJFNVTLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
