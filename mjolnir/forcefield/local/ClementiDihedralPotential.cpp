#include <mjolnir/forcefield/local/ClementiDihedralPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ClementiDihedralPotential<double>;
template class ClementiDihedralPotential<float>;
} // mjolnir
