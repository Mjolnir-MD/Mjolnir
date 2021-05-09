#include <mjolnir/omp/GFWNpTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
