#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ProteinDNANonSpecificInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ProteinDNANonSpecificInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ProteinDNANonSpecificInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ProteinDNANonSpecificInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
