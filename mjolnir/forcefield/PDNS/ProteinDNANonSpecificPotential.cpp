#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ProteinDNANonSpecificPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ProteinDNANonSpecificPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ProteinDNANonSpecificPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ProteinDNANonSpecificPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
