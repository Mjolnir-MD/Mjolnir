#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ProteinDNANonSpecificParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ProteinDNANonSpecificParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ProteinDNANonSpecificParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ProteinDNANonSpecificParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
