#include <mjolnir/core/GFWNpTLangevinIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GFWNpTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class GFWNpTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class GFWNpTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class GFWNpTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
