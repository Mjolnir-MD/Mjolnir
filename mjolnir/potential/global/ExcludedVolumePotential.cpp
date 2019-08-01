#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ExcludedVolumePotential<double>;
template class ExcludedVolumePotential<float>;
} // mjolnir
