#include <mjolnir/forcefield/FLP/FlexibleLocalAnglePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class FlexibleLocalAnglePotential<double>;
template class FlexibleLocalAnglePotential<float>;
} // mjolnir
