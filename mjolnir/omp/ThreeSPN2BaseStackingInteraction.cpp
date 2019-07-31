#include <mjolnir/omp/ThreeSPN2BaseStackingInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ThreeSPN2BaseStackingInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BaseStackingInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BaseStackingInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BaseStackingInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
