#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteractionPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2BaseBaseInteractionPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteractionPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteractionPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BaseBaseInteractionPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
