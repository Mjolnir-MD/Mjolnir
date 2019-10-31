#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
