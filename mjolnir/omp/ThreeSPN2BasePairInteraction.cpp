#include <mjolnir/omp/ThreeSPN2BasePairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ThreeSPN2BasePairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BasePairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BasePairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BasePairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
