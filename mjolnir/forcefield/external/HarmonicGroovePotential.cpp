#include <mjolnir/forcefield/external/HarmonicGroovePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class HarmonicGroovePotential<double>;
template class HarmonicGroovePotential<float>;
} // mjolnir
