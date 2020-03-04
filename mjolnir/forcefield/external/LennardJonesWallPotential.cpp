#include <mjolnir/forcefield/external/LennardJonesWallPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LennardJonesWallPotential<double>;
template class LennardJonesWallPotential<float>;
} // mjolnir
