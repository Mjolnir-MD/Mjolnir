#include <mjolnir/omp/ThreeSPN2BaseBaseInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
