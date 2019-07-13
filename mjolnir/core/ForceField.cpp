#include <mjolnir/core/ForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ForceField<SimulatorTraits<double, UnlimitedBoundary>>;
template class ForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
template class ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
