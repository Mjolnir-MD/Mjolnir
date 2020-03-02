#include <mjolnir/potential/local/HarmonicPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class HarmonicPotential<double>;
template class HarmonicPotential<float>;
} // mjolnir
