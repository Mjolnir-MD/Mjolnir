#include <mjolnir/omp/gBAOABLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class gBAOABLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class gBAOABLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class gBAOABLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class gBAOABLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
