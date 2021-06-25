#include <mjolnir/forcefield/global/LennardJonesAttractivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LennardJonesAttractivePotential<double>;
template class LennardJonesAttractivePotential<float >;
} // mjolnir
