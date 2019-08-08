#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
