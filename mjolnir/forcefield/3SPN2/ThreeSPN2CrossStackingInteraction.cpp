#include <mjolnir/forcefield/3SPN2/ThreeSPN2CrossStackingInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
