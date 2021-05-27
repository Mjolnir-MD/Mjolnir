#include <mjolnir/omp/GFWNPTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GFWNPTLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class GFWNPTLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class GFWNPTLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GFWNPTLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
