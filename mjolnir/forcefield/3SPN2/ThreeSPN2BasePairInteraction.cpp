#include <mjolnir/forcefield/3SPN2/ThreeSPN2BasePairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ThreeSPN2BasePairInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BasePairInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BasePairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BasePairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
