#include <mjolnir/omp/BAOABLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
