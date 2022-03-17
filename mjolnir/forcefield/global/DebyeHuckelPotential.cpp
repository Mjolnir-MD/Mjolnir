#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class DebyeHuckelParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class DebyeHuckelParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class DebyeHuckelParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class DebyeHuckelParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
