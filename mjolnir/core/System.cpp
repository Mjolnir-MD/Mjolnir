#include <mjolnir/core/System.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class System<SimulatorTraits<double, UnlimitedBoundary>>;
template class System<SimulatorTraits<float,  UnlimitedBoundary>>;
template class System<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class System<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
