#include <mjolnir/omp/ThreeSPN2CrossStackingInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2CrossStackingInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2CrossStackingInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2CrossStackingInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2CrossStackingInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
