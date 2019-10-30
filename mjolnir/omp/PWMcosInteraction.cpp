#include <mjolnir/omp/PWMcosInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class PWMcosInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class PWMcosInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
template class PWMcosInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class PWMcosInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
