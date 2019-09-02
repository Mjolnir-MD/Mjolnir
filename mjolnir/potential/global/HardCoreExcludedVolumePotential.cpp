#include <mjolnir/potential/global/HardCoreExcludedVolumePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class HardCoreExcludedVolumePotential<double>;
template class HardCoreExcludedVolumePotential<float >;
} // mjolnir
