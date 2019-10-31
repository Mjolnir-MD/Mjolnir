#include <mjolnir/omp/ProteinDNANonSpecificInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ProteinDNANonSpecificInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class ProteinDNANonSpecificInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
template class ProteinDNANonSpecificInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ProteinDNANonSpecificInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
