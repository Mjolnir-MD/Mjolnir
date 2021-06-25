#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteractionPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2BaseBaseInteractionParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteractionParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteractionParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BaseBaseInteractionParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
