#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseStackingInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
