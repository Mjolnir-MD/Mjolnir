#include <mjolnir/core/GFWNPTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GFWNPTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class GFWNPTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class GFWNPTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GFWNPTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
