#include <mjolnir/core/VelocityVerletIntegrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
template class VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
template class VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
