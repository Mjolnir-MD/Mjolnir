#include <mjolnir/forcefield/AFMFit/AFMFitInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class AFMFitInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class AFMFitInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class AFMFitInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class AFMFitInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
