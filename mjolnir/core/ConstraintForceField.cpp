#include <mjolnir/core/ConstraintForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ConstraintForceField<SimulatorTraits<double, UnlimitedBoundary>>;
template class ConstraintForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
template class ConstraintForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ConstraintForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
