#include <mjolnir/potential/external/ExcludedVolumeWallPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ExcludedVolumeWallPotential<double>;
template class ExcludedVolumeWallPotential<float>;
} // mjolnir
