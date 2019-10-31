#include <mjolnir/forcefield/PWMcos/PWMcosInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class PWMcosInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class PWMcosInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class PWMcosInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class PWMcosInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
