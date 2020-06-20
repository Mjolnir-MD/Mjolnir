#include <mjolnir/input/read_constraint.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ConstraintForceField<SimulatorTraits<double, UnlimitedBoundary>>
    read_constraint(const toml::value& constraint);
template ConstraintForceField<SimulatorTraits<float,  UnlimitedBoundary>>
    read_constraint(const toml::value& constarint);
template ConstraintForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>
    read_constraint(const toml::value& constraint);
template ConstraintForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>
    read_constraint(const toml::value& constraint);
} // mjolnir
