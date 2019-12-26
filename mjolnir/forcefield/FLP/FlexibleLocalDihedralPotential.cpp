#include <mjolnir/forcefield/FLP/FlexibleLocalDihedralPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class FlexibleLocalDihedralPotential<double>;
template class FlexibleLocalDihedralPotential<float>;
} // mjolnir
