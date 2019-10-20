#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ProteinDNANonSpecificPotential<double>;
template class ProteinDNANonSpecificPotential<float >;

} // mjolnir
